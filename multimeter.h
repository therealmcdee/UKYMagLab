#pragma once

#include <visa.h>

class MultiMeter
{
private:
	const char zVISAAddress[43] = "USB0::0x0957::0x0A07::MY48001873::0::INSTR";
	const char yVISAAddress[43] = "USB0::0x0957::0x0A07::MY48001709::0::INSTR";
	const char xVISAAddress[43] = "USB0::0x0957::0x0A07::MY48002121::0::INSTR";
    const char pwrVISAAddress[13] = "ASRL5::INSTR";  // power supply for the inner coil 
    const char pwrVISAAddress2[13] = "ASRL3::INSTR"; // power supply for outer coil 
    const char pwrVISAAddress3[13] = "ASRL4::INSTR"; // third power supply for the gradient coil 

    ViConstRsrc yDMM = xVISAAddress;
    ViConstRsrc xDMM = yVISAAddress;
    ViConstRsrc zDMM = zVISAAddress;
    ViConstRsrc pwr  = pwrVISAAddress;
    ViConstRsrc pwr2  = pwrVISAAddress2;
    ViConstRsrc pwr3  = pwrVISAAddress3; // GC

    ViStatus Status;

    ViSession defaultRMx, vix;
    ViSession defaultRMy, viy;
    ViSession defaultRMz, viz;
    ViSession defaultRMpwr, vipwr;
    ViSession defaultRMpwr2, vipwr2;
    ViSession defaultRMpwr3, vipwr3; //GC

public:
    MultiMeter()
    {
        Status = viOpenDefaultRM(&defaultRMx);
        Status |= viOpen(defaultRMx, xDMM, VI_NULL, VI_NULL, &vix);

        Status = viOpenDefaultRM(&defaultRMy);
        Status |= viOpen(defaultRMy, yDMM, VI_NULL, VI_NULL, &viy);

        Status = viOpenDefaultRM(&defaultRMz);
        Status |= viOpen(defaultRMz, zDMM, VI_NULL, VI_NULL, &viz);

        Status = viOpenDefaultRM(&defaultRMpwr);
        Status |= viOpen(defaultRMpwr, pwr, VI_NULL, VI_NULL, &vipwr);

        Status = viOpenDefaultRM(&defaultRMpwr2);
        Status |= viOpen(defaultRMpwr2, pwr2, VI_NULL, VI_NULL, &vipwr2);

        Status = viOpenDefaultRM(&defaultRMpwr3);                               //GC
        Status |= viOpen(defaultRMpwr3, pwr3, VI_NULL, VI_NULL, &vipwr3);       //GC


        Status |= viPrintf(vix, "*RST\n");
        Status |= viPrintf(viy, "*RST\n");
        Status |= viPrintf(viz, "*RST\n");



        Status |= viPrintf(vipwr, "VOLT 50\n");
        Status |= viPrintf(vipwr2, "VOLT 50\n");
        Status |= viPrintf(vipwr3, "VOLT 35\n");                              //GC

        //  change current 

       
        setCurrent(1.25, 1.25, 1.25); // change the current ratio between outer and inner coils 
 
  
        //Status |= viPrintf(vipwr, "CURR 1\n");
        //Status |= viPrintf(vipwr, "APPL:SOUR:CURR 1\n");

        Status |= viPrintf(vix, "SENS:VOLT:DC:NPLC 2\n");
        Status |= viPrintf(viy, "SENS:VOLT:DC:NPLC 2\n");
        Status |= viPrintf(viz, "SENS:VOLT:DC:NPLC 2\n");

        Status |= viPrintf(vix, "VOLT:DC:IMP:AUTO 1\n");
        Status |= viPrintf(viy, "VOLT:DC:IMP:AUTO 1\n");
        Status |= viPrintf(viz, "VOLT:DC:IMP:AUTO 1\n");

        Status |= viPrintf(vix, "VOLT:RANG 10\n");
        Status |= viPrintf(viy, "VOLT:RANG 10\n");
        Status |= viPrintf(viz, "VOLT:RANG 10\n");

        Status |= viPrintf(vix, "TRIG:SOUR IMM\n");
        Status |= viPrintf(viy, "TRIG:SOUR IMM\n");
        Status |= viPrintf(viz, "TRIG:SOUR IMM\n");

    }

    /**
    void setCurrent(float curr,float curr2) // when 2 powerSpply
    {
        Status |= viPrintf(vipwr, "CURR %.4f\n", curr);
        Status |= viPrintf(vipwr2, "CURR %.4f\n", curr2);
    } **/

 
    void setCurrent(float curr, float curr2, float curr3) // 3 PSupply
    {
        Status |= viPrintf(vipwr, "CURR %.4f\n", curr);
        Status |= viPrintf(vipwr2, "CURR %.4f\n", curr2);
        Status |= viPrintf(vipwr3, "CURR %.4f\n", curr3);
    }  


    void outputOn(const int& state)
    {
        if (state == 1)
        {
            Status |= viPrintf(vipwr, "OUTP ON\n");
            Status |= viPrintf(vipwr2, "OUTP ON\n");
            Status |= viPrintf(vipwr3, "OUTP ON\n"); // GC
        }
        else if (state == 0)
        {
            Status |= viPrintf(vipwr, "OUTP OFF\n");
            Status |= viPrintf(vipwr2, "OUTP OFF\n");
            Status |= viPrintf(vipwr3, "OUTP OFF\n"); //GC
        }

        if (Status != 0)
            std::cout << "outputOn : I/O Error" << std::endl;
    }

    double measureBx()
    {
        double res = 0;

        for(int i = 0; i < 5; i++)
        {
            double temp;
            Status |= viPrintf(vix, "READ?\n");
            Status |= viScanf(vix, "%lf", &temp);
            res += temp / 5.0;
        }
        return 100. * res;
    }



    double measureBy()
    {
        double res = 0;

        for(int i = 0; i < 5; i++)
        {
            double temp;
            Status |= viPrintf(viy, "READ?\n");
            Status |= viScanf(viy, "%lf", &temp);
            res += temp / 5.0;
        }
        return 100. * res;
    }

    double measureBz()
    {
        double res = 0;

        for(int i = 0; i < 5; i++)
        {
            double temp;
            Status |= viPrintf(viz, "READ?\n");
            Status |= viScanf(viz, "%lf", &temp);
            res += temp / 5.0;
        }
        return 100. * res;
    }

	void getVISA_Addr()
	{
		std::cout << " X : " << xVISAAddress << std::endl;
		std::cout << " Y : " << yVISAAddress << std::endl;
		std::cout << " Z : " << zVISAAddress << std::endl;
	}

    ~MultiMeter()
    {
        Status |= viClose(vix);
        Status |= viClose(defaultRMx);
        Status |= viClose(viy);
        Status |= viClose(defaultRMy);
        Status |= viClose(viz);
        Status |= viClose(defaultRMz);
        Status |= viClose(vipwr);
        Status |= viClose(defaultRMpwr);

        Status |= viClose(vipwr2);
        Status |= viClose(defaultRMpwr2);

        Status |= viClose(vipwr3);
        Status |= viClose(defaultRMpwr3);
    }
};
