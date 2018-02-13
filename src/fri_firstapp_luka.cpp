/*----------------------------------------------------------------
KUKA control
HRI2 lab, IIT
----------------------------------------------------------------*/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <math.h>
#include <limits.h>
#include <friudp.h>
#include <friremote.h>
#include <armadillo>


#include <signal.h>
#include <sys/mman.h>
#include <execinfo.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef __XENO__
    #include <rtdk.h>
#endif

#include <boost/circular_buffer.hpp>
//#include <boost/tokenizer.hpp>

#include <ati_iface.h>

#define NSEC_PER_SEC	1000000000ULL

#define EPSilon 2.2204460492503131e-016

inline uint64_t get_time_ns(clockid_t clock_id=CLOCK_MONOTONIC)
{
    uint64_t time_ns;
    struct timespec ts;
    clock_gettime(clock_id, &ts);
    time_ns = ts.tv_sec * Nclude <stdlirtdkb.h>
SEC_PER_SEC + ts.tv_nsec;
    return time_ns;
}

static int run_loop = 1;

typedef struct {
    uint64_t    ts;
    float       ati[6];
    void sprint(char *buff, size_t size) {
        snprintf(buff, size, "%lu\t%f\t%f\t%f\t%f\t%f\t%f\t0\t0\t0\t0\t0\t0\n", ts,
                 ati[0],ati[1],ati[2],ati[3],ati[4],ati[5]);
    }
    void fprint(FILE *fp) {
        fprintf(fp, "%lu\t%f\t%f\t%f\t%f\t%f\t%f\t0\t0\t0\t0\t0\t0\n", ts,
                 ati[0],ati[1],ati[2],ati[3],ati[4],ati[5]);
    }
} sens_data_t ; // FT COMMENT


static void shutdown(int sig __attribute__((unused)))
{
    run_loop = 0;
    printf("got signal .... Shutdown\n");
}

static void set_signal_handler(void)
{
    signal(SIGINT, shutdown);
    signal(SIGINT, shutdown);
    signal(SIGKILL, shutdown);
}

// *** ***********************

#ifndef M_PI
#define M_PI 3.14159
#endif

using namespace std;
using namespace arma;

bool runOnce = true, nullControlFlag = false, startFlag = true;
float deltaT = 0.005; //0.02;


#include <stdio.h>
#include <sys/select.h>
#include <termios.h>
#include <stropts.h>
#include <sys/ioctl.h>

int _kbhit() {
    static const int STDIN = 0;
    static bool initialized = false;

    if (! initialized) {
        // Use termios to turn off line buffering
        termios term;
        tcgetattr(STDIN, &term);
        term.c_lflag &= ~ICANON;
        tcsetattr(STDIN, TCSANOW, &term);
        setbuf(stdin, NULL);
        initialized = true;
    }

    int bytesWaiting;
    ioctl(STDIN, FIONREAD, &bytesWaiting);
    return bytesWaiting;
}

#define DATASIZE1        1e7
int trialnumbertext=11;

char textfilename[50];

int
#ifdef __WIN32

_tmain
#else
#ifdef _WRS_KERNEL
friFirstApp
#else
main
#endif
#endif

(int argc, char *argv[])
{
    // sampling
    float *Positionar = (float*) malloc(DATASIZE1*sizeof(float));
    int myPocounter = 0;


    char key;


    // ***** /////////////////////////////////////
    boost::circular_buffer<sens_data_t> sens_log;
    sens_log.set_capacity(LOG_SIZE);

    set_signal_handler();  // FT COMMENT

    ///////////////////////////////////////////////////////////////////////
    Ati_Sens * ati = new Ati_Sens(true); // FT COMMENT
    ///////////////////////////////////////////////////////////////////////

    ati_log_t   sample;
    sens_data_t sens_data, sens_data_bias; // FT COMMENT
    uint64_t    start = get_time_ns();
    // ***** /////////////////////////////////////



    cout << "Opening FRI Interface" << endl;
    {

        friRemote friInst(49948,"192.168.0.10" );
        //friRemote friInst(49938,"192.168.0.20" );
        FRI_QUALITY lastQuality = FRI_QUALITY_BAD;
        FRI_CTRL lastCtrlScheme = FRI_CTRL_OTHER;

        float newJntVals_msrd[LBR_MNJ];
        colvec measured_for(LBR_MNJ,1);


        float measCartVals[FRI_CART_FRM_DIM];
        float newCartVals[FRI_CART_FRM_DIM];
        float newForceTorqueAdd[FRI_CART_VEC];
        float myCartstiffness[6]={1500.0,1500.0,1500.0,150.0,150.0,150.0};
        float myJntstiffness[7]={0,0,0,0,0,0,0};






        printf("START\n");
        friInst.doDataExchange();

        while(!(key = _kbhit()))
        {

            friInst.setToKRLInt(0,1);
            if ( friInst.getQuality() >= FRI_QUALITY_OK)
            {
                friInst.setToKRLInt(0,10);
            }

            friInst.setToKRLReal(0,friInst.getFrmKRLReal(1));

            if ( lastCtrlScheme != friInst.getCurrentControlScheme())
            {
                cout << "switching control scheme " << lastCtrlScheme;
                lastCtrlScheme = friInst.getCurrentControlScheme();
                cout << " to " << lastCtrlScheme<<endl;
            }

            // FORCE SENSOR
            ati->get_last_sample(sample);

            sens_data.ts = get_time_ns() - start;
            memcpy((void*)&sens_data.ati, &sample.ft, sizeof(float)*6);
            sens_log.push_back(sens_data);
            //sens_data.fprint(stderr);

            measured_for<<(sens_data.ati[0] - sens_data_bias.ati[0])/1000<<(sens_data.ati[1] - sens_data_bias.ati[1])/1000<<(sens_data.ati[2] - sens_data_bias.ati[2])/1000<<(sens_data.ati[3] - sens_data_bias.ati[3])/1000<<(sens_data.ati[4] - sens_data_bias.ati[4])/1000<<(sens_data.ati[5] - sens_data_bias.ati[5])/1000;
            // FT COMMENT

            //printf("F/T:\n");
            //measured_for.print();

            switch (friInst.getCurrentControlScheme())
            {
            case   FRI_CTRL_POSITION:
            case   FRI_CTRL_JNT_IMP:
                {

                    // CONTROL

                    for (int i = 0; i < LBR_MNJ; i++) // measured joint values
                    {
                        newJntVals_msrd[i] = friInst.getMsrMsrJntPosition()[i];
                    }

                    for (int i = 0; i < FRI_CART_FRM_DIM; i++) // measured Cartesian values
                    {
                        measCartVals[i] = friInst.getMsrCartPosition()[i];
                        newCartVals[i] = friInst.getMsrCmdCartPosition()[i];
                    }

                    //cout << measCartVals[3] << "\t" << measCartVals[7] << "\t" << measCartVals[11] << " " << endl;
                    //cout << measured_for[0] << "\t" << measured_for[1] << "\t" << measured_for[2] << " " << endl;


                    if ( friInst.getState() == FRI_STATE_CMD)
                    {
                        if ( friInst.isPowerOn() )
                        {

                            // impedance controller here


                            Positionar[myPocounter]=measCartVals[0]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[1]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[2]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[3]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[4]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[5]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[6]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[7]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[8]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[9]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[10]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measCartVals[11]; // z
                            myPocounter++;
                            Positionar[myPocounter]=measured_for[0]; // ft
                            myPocounter++;
                            Positionar[myPocounter]=measured_for[1]; // ft
                            myPocounter++;
                            Positionar[myPocounter]=measured_for[2]; // ft
                            myPocounter++;
                            Positionar[myPocounter]=measured_for[3]; // ft
                            myPocounter++;
                            Positionar[myPocounter]=measured_for[4]; // ft
                            myPocounter++;
                            Positionar[myPocounter]=measured_for[5]; // ft
                            myPocounter++;

                        }

                    }

                    friInst.doJntImpedanceControl(NULL,myJntstiffness,NULL,NULL); //newJntAddTorque

                }
                break;

            case FRI_CTRL_CART_IMP:

                    {

                        for (int i = 0; i < LBR_MNJ; i++) // measured joint values
                        {
                            newJntVals_msrd[i] = friInst.getMsrMsrJntPosition()[i];
                        }

                        for (int i = 0; i < FRI_CART_FRM_DIM; i++) // measured Cartesian values
                        {
                            measCartVals[i] = friInst.getMsrCartPosition()[i];
                            newCartVals[i] = friInst.getMsrCmdCartPosition()[i];
                        }



                        if ( friInst.getState() == FRI_STATE_CMD)
                        {
                            if ( friInst.isPowerOn() )
                            {

                                // CARTESIAN POSITION
                                //newCartVals[3] = newCartVals[3];
                                //newCartVals[7] = newCartVals[7];
                                //newCartVals[11] = newCartVals[11];


                                Positionar[myPocounter]=measCartVals[0]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[1]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[2]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[3]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[4]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[5]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[6]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[7]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[8]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[9]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[10]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measCartVals[11]; // z
                                myPocounter++;
                                Positionar[myPocounter]=measured_for[0]; // ft
                                myPocounter++;
                                Positionar[myPocounter]=measured_for[1]; // ft
                                myPocounter++;
                                Positionar[myPocounter]=measured_for[2]; // ft
                                myPocounter++;
                                Positionar[myPocounter]=measured_for[3]; // ft
                                myPocounter++;
                                Positionar[myPocounter]=measured_for[4]; // ft
                                myPocounter++;
                                Positionar[myPocounter]=measured_for[5]; // ft
                                myPocounter++;

                            }

                        }

                        //friInst.doCartesianImpedanceControl(newCartVals,myCartstiffness,NULL,NULL);
                        friInst.doCartesianImpedanceControl(NULL,NULL,NULL,NULL);
                    }
                    break;

            // Cartesian impedance control case was here

            default:
                friInst.doDataExchange();
            }

            if ( friInst.getFrmKRLInt(0) == -1)
            {
                cout << "leaving \n";
                break;
            }

            if ( friInst.getQuality() != lastQuality)
            {
                cout << "quality change detected "<< friInst.getQuality()<< " \n";
                cout << friInst.getMsrBuf().intf;
                cout << endl;
                lastQuality=friInst.getQuality();
            }
        }

    }

    // SAVE DATA TO FILE
    int wctr=0;
    ofstream arrayofPo;
    sprintf(textfilename,"trial%d.txt",trialnumbertext);

    arrayofPo.open (textfilename);

    while (wctr<DATASIZE1)
    {
       if (abs(Positionar[wctr])>1e7) break;

       vec test_norm;
       test_norm<<Positionar[wctr]<<Positionar[wctr+1]<<Positionar[wctr+2]<<Positionar[wctr+3]<<Positionar[wctr+4]<<Positionar[wctr+5]<<Positionar[wctr+6]<<Positionar[wctr+7]<<Positionar[wctr+8]<<Positionar[wctr+9]<<Positionar[wctr+10]<<Positionar[wctr+11]<<Positionar[wctr+12]<<Positionar[wctr+13]<<Positionar[wctr+14]<<Positionar[wctr+15]<<Positionar[wctr+16]<<Positionar[wctr+17];
       if (norm(test_norm) == 0) break;

       arrayofPo<<Positionar[wctr]<<" "<<Positionar[wctr+1]<<" "<<Positionar[wctr+2]<<" "<<Positionar[wctr+3]<<" "<<Positionar[wctr+4]<<" "<<Positionar[wctr+5]<<" "<<Positionar[wctr+6]<<" "<<Positionar[wctr+7]<<" "<<Positionar[wctr+8]<<" "<<Positionar[wctr+9]<<" "<<Positionar[wctr+10]<<" "<<Positionar[wctr+11]<<" "<<Positionar[wctr+12]<<" "<<Positionar[wctr+13]<<" "<<Positionar[wctr+14]<<" "<<Positionar[wctr+15]<<" "<<Positionar[wctr+16]<<" "<<Positionar[wctr+17]<<endl;
       wctr=wctr + 18;
    }
    arrayofPo.close();
    delete(Positionar);



    return EXIT_SUCCESS;
}

/*

bool Cmove(const mat &T, Real duration) {
    cout<<"***Start Cmove***\n";
    SetStrategy(2);

    std::ofstream fout("CmoveT.txt");
    //!     double dt = 1.0/SERVO_RATE;
    int N = int(duration/dt);
    Transform3D Tb = World2Base(T,(*this).WhichArm);
    //!starts from the last commanded position
    Transform3D T0 = getRobotT('c','b');
    Transform3D *Ti;
    Ti = Transform3D::ctraj(T0,Tb,duration/dt);

    //! Execute
    run();
    for (int i=0; i<N; i++){
        setCartPos(Ti[i]);
        send();
        fout << QFRAME(Ti[i]);
    }
    stop();
    fout.flush();   fout.close();

    //!update commanded joint position
    updateCmdJnt();
    cout<<"\nDone CmoveT...\n"<<"***************\n";

    return true;
}


mat * ctraj(const mat & T0, const mat & T1, double n)
     {
        if ( (int)n <= 0)
            cout<<"\nInvalid n for ctraj\n";
        mat *Ti = new mat[(int)n+1];
        colvec i((int)n), r((int)n);
        for (int k=1; k<=(int)n; k++){
            Ti[k-1] = mat(4,4);
            Ti[k-1](4,4) = 0.0;
            i(k) = static_cast<double>(k);
            if (k>1)
                r(k) = double(i(k-1)/(n-1));
        }
        for (int k=0; k<(int)n; k++){
            Ti[k] = trinterp(T0, T1, r(k+1));
        }
        Ti[(int)n] = T1;
        return Ti;
    }

mat trinterp(const mat & T0, const mat & T1, double r)
    {
        colvec p0(3,1), p1(3,1), pr(3,1);
        colvec q0, q1, qr;
        mat Rr;

        q0 = R2Q(T0);
        q1 = R2Q(T1);

        p0 = T0.P();
        p1 = T1.P();

        qr = qinterp(q0, q1, r);
        pr = p0*(1-r) + r*p1;

        Rr = Quaternion::Q2R(qr);
        return Transform3D(pr, Rr);
    }


colvec qinterp(colvec Q1, colvec Q2, double r)
    {
        if ((r < 0) || (r > 1))
            cerr << "qinterp(Q1, Q2, r): r < 0 or r > 1. r is set to 0." << endl;
        Real theta = Q1(0)*Q2(0) + Q1(1)*Q2(1) + Q1(2)*Q2(2) + Q1(3)*Q2(3);
        theta = acos( (theta > 1 ? 1 : (theta < -1 ? -1 : theta)) );
        colvec q;
        if (abs(theta) < EPSilon)
            q = Q1;
        else
            q = (sin((1-r)*theta) * Q1 + sin(r*theta) * Q2) / sin(theta) ;
        return q;
    }
colvec R2Q(const mat &R)
     {
        colvec Q;
        double s,x,y,z;
        double ss,n,x1,y1,z1;
        bool c;
        s = R(0,0) + R(1,1) + R(2,2) + 1.0;
        if (s < 0.0)
            s = 0.0;
        else
            s = static_cast<double>(0.5) * sqrt( s );
        x = R(3-1,2-1) - R(2-1,3-1);
        y = R(1-1,3-1) - R(3-1,1-1);
        z = R(2-1,1-1) - R(1-1,2-1);
        if ( (R(0,0) >= R(1,1)) & (R(0,0) >= R(2,2)) ){
            x1 = R(0,0) - R(1,1) - R(2,2) + 1.0;
            y1 = R(2-1,1-1) + R(1-1,2-1);
            z1 = R(3-1,1-1) + R(1-1,3-1);
            c = (x >= 0);
        }else if ( R(2-1,2-1) >= R(3-1,3-1) ){
            x1 = R(2-1,1-1) + R(1-1,2-1);
            y1 = R(2-1,2-1) - R(1-1,1-1) - R(3-1,3-1) + 1.0;
            z1 = R(3-1,2-1) + R(2-1,3-1);
            c = (y >= 0);
        }else{
            x1 = R(3-1,1-1) + R(1-1,3-1);
            y1 = R(3-1,2-1) + R(2-1,3-1);
            z1 = R(3-1,3-1) - R(1-1,1-1) - R(2,2) + ONE;
            c = (z >= 0);
        }
        if (c){
            x += x1;
            y += y1;
            z += z1;
        }else{
            x -= x1;
            y -= y1;
            z -= z1;
        }
        n = sqrt( x*x + y*y + z*z );
        if (n == 0){
            Q(0) + 1.0;
            Q(1) = 0;       Q(2) = 0;       Q(3) = 0;
         }else{
            colvec v(3);
            Q(0) + s;
            ss = sqrt( 1 - s*s )/n;
            Q(1)=ss*x;      Q(2)=ss*y;      Q(3)=ss*z;
        }
        return Q;
    }
    */
