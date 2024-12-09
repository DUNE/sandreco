#ifndef TigerElectronics_h
#define TigerElectronics_h
#include "Common.h"
#include "ElectronicChannel.h"
using namespace std;
namespace TIGER {
    class Readout {
    public:
        //Constructor
        Readout();  //int setup, int electronics, Geometry* geometry
        //Destructor
        ~Readout();
        //Function
        bool                       Get_PrintPng() { return PrintPng; };
        bool                       Get_PrintPDF() { return PrintPDF; };
        void                       Simulate_TIGER(); //vector<double>
        void                       Induce_on_channel(vector<double> wf);
        void                       Integrate_Charge();
        void                       Set_ChannelID(int io) { channel_id = io; };
        void                       Set_PrintPng(bool io) { PrintPng = io; };
        void                       Set_PrintPDF(bool io) { PrintPDF = io; };
        ElectronicChannel*         Get_Channel() { return channel; };
        vector<double>             GetDigitalOutput() { return output; };

    private:
        //Variable
        //bool      PrintInfo;
        bool      PrintPDF;
        bool      PrintPng;
        int       channel_id;
        ElectronicChannel* channel;  
        vector<double> output;

        //Function
        //int       Get_ChannelID(int type);
        void      Initialize_TIGER();
        void      Integration_TIGER();
        void      Extract_Charge_Time();
        double    Get_Charge_TIGER(ElectronicChannel* ch);
        double    Get_Time_TIGER(ElectronicChannel* ch);
        double    Get_dTime_TIGER(ElectronicChannel* ch);

        TF1* f[n_ns];
        //APV
        double time, dtime;
        void Set_Time(double io) { time = io; };
        double Get_Time() { return time; };
        void Set_dTime(double io) { dtime = io; };
        double Get_dTime() { return dtime; };
        //Tiger
        double T_branch(double t) {
            t = t - 10.;
            if (t > 0) return 2000 * (0.00181928 * exp(-t / 3.) - 0.0147059 * exp(-t / 20.) + 0.0128866 * exp(-t / 100.));
            return 0;
        };
        double E_branch(double t) {
            t = t - 5;
            if (t > 0) return 0.000627357 * (1358.7 * exp(-t * 0.0385647) * t + 1358.7 * exp(-t * 0.0114353) * t + 100164. * exp(-t * 0.0385647) - 100164. * exp(-0.0114353 * t));
            return 0;
        };

    };
}
#endif