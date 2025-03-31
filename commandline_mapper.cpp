#include <iostream>
#include <string>
#include "filereader.h"
extern "C" {
#include "NIDAQmx.h"
}
#define _USE_MATH_DEFINES
#include <visa.h>
//#include "multimeter.h"
#include "pubSysCls.h"
#include <fstream>
#include <iomanip>
#include "Eigen/Eigen"
#include "Eigen/QR"
#include <math.h>
#include <vector>
#include <Windows.h>


using Eigen::MatrixXd;
using Eigen::VectorXd;

bool isSubtractionOn = false;




// Send message and wait for newline
void msgUser(const char* msg) {
	std::cout << msg;
	getchar();
}

float findBestFreq(float64* data, float w_init) {
	float w_best = 0;
	int NumSamples = 5000;
	int SampleRate = 5000;

	int numVals = 5;
	MatrixXd w = MatrixXd::Zero(NumSamples, 3);
	VectorXd t = VectorXd::LinSpaced(NumSamples, 0, NumSamples / SampleRate);
	VectorXd wVals = VectorXd::LinSpaced(numVals, 0.99f * w_init, 1.01f * w_init); 
	VectorXd chisq = VectorXd::Zero(numVals);
	
	float w_test = 0;
	for (int i = 0; i < numVals; i++) {
		w_test = wVals(i);
		w.col(0) = (2 * M_PI * t * w_test).array().sin();
		w.col(1) = (2 * M_PI * t * w_test).array().cos();
		w.col(2) = VectorXd::Ones(NumSamples);
		MatrixXd wt = w.transpose();
		MatrixXd pinv = (wt * w).inverse() * wt;

		
		Eigen::Map<MatrixXd>   Bt(data, 4, NumSamples);
		MatrixXd B = Bt.transpose();
		MatrixXd out = pinv * B;
		//VectorXd v = w * out.col(3);
		VectorXd v = B.col(3) - w * out.col(3);
		chisq(i) = v.transpose() * v;
		
	}
	
	MatrixXd w2 = MatrixXd::Zero(numVals, 3);
	w2.col(0) = wVals.array().square();
	w2.col(1) = wVals;
	w2.col(2) = VectorXd::Ones(numVals);

	MatrixXd wt2 = w2.transpose();
	MatrixXd pinv2 = (wt2 * w2).inverse() * wt2;

	MatrixXd out2 = pinv2 * chisq;

	/**/
	
	std::cout << "--------------------------" << std::endl;
	//std::cout << wVals << std::endl;
	std::cout << chisq << std::endl;
	std::cout << "Best Freq: " << -1*out2(1) / (2 * out2(0)) << std::endl;
	std::cout << "--------------------------" << std::endl;

	w_best = -1 * out2(1) / (2 * out2(0));
	return w_best;
}

ViStatus powerOn(ViSession& pwr, ViStatus Status, bool on) {
	if (on) {
		Status |= viPrintf(pwr, "OUTP ON\n");
	}
	else {
		Status |= viPrintf(pwr, "OUTP OFF\n");
	}
	return Status;
}

//ViStatus powerOn(ViSession& pwr, ViStatus Status) {
//	Status |= viPrintf(pwr, "OUTP ON\n");
//	return Status;
//}

void setCurrent(ViSession& pwr, ViStatus Status, float curr) {
	Status |= viPrintf(pwr, "CURR %.4f\n", curr);
}

void printVec3dList(const std::vector<Vec3d>& list)
{
	for (const Vec3d& p : list)
		std::cout << p << std::endl;
}

/*
void measureAndPrint(MultiMeter& mm)
{
	std::cout << mm.measureBx() << " " << mm.measureBy() << " " << mm.measureBz() << std::endl;
}
*/

/*
void togglePWR(MultiMeter& mm, int s)
{
	mm.outputOn(s);
}
*/

Vec3d measureField(TaskHandle& taskHandle, int NumSamples, int32* read, float64* data) {

	DAQmxReadAnalogF64(taskHandle, NumSamples, 10.0, DAQmx_Val_GroupByScanNumber, data, 4 * NumSamples, read, NULL);
	// For now, let's average over all of the samples we just collected to return a single Vec3d measurement
	float64 Bx = 0, By = 0, Bz = 0, V = 0;
	for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
		//write [t,Bx,By,Bz] to the file for each line
		Bx += data[i] / NumSamples;
		By += data[i + 1] / NumSamples;
		Bz += data[i + 2] / NumSamples;
		V += data[i + 3] / NumSamples;
	}
	Vec3d measurementValue = {Bx,By,Bz};

	return measurementValue;
}

void testMeasureField(TaskHandle& taskHandle, int NumSamples) {

	int32 read = 0;
	// initialize the array to hold the x,y,z data
	float64* data = new float64[4 * NumSamples];
	//Sleep(5 * 1e3);
	DAQmxReadAnalogF64(taskHandle, NumSamples, 10.0, DAQmx_Val_GroupByScanNumber, data, 4 * NumSamples, &read, NULL);
	// For now, let's average over all of the samples we just collected to return a single Vec3d measurement
	float64 Bx = 0, By = 0, Bz = 0, V = 0;
	for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
		//write [t,Bx,By,Bz] to the file for each line
		Bx += data[i] / NumSamples;
		By += data[i + 1] / NumSamples;
		Bz += data[i + 2] / NumSamples;
		V += data[i + 3] / NumSamples;
	}

	std::cout << "Bx: " << Bx 
		<< ", By: " << By 
		<< ", Bz: " << Bz
		<< ", I: " << V 
		<< std::endl;
	
	/*
	DAQmxReadAnalogF64(taskHandle, NumSamples, 10.0, DAQmx_Val_GroupByScanNumber, data, 3 * NumSamples, &read, NULL);
	// For now, let's average over all of the samples we just collected to return a single Vec3d measurement
	float64 Bx = 0, By = 0, Bz = 0;
	for (uInt32 i = 0; i < 3 * NumSamples; i += 3) {
		//write [t,Bx,By,Bz] to the file for each line
		Bx += data[i] / NumSamples;
		By += data[i + 1] / NumSamples;
		Bz += data[i + 2] / NumSamples;
	}
	std::cout << "Bx: " << Bx << ", By: " << By << ", Bz: " << Bz << "\n";
	*/
	delete[] data;
	
}



int Map(const std::vector<Vec3d>& list, TaskHandle taskHandle, int NumSamples, std::vector<std::vector<double>>& result, const Eigen::Ref<const Eigen::MatrixXd>& pinv3, int back_flag, ViSession& vipwr, ViStatus Status)
{
	//NIDAQmx stuff
	int32 read = 0;
	// initialize the array to hold the x,y,z data
	float64* data = new float64[4 * NumSamples];
	
	//std::vector<Vec3d> emptyResult;
	std::vector<std::vector<double>> emptyResult;

	result = emptyResult;
	Vec3d probePosition{ 0, 0, 0 };

	int xyCmToCounts = 750;
	int zCmToCounts = 3000;

	int ACC_LIM_RPM_PER_SEC = 30;
	int VEL_LIM_RPM = 50;
	int TIME_TILL_TIMEOUT = 10000;
	int MOVE_DISTANCE_CNTS;

	size_t portCount = 0;
	std::vector<std::string> comHubPorts;

	//Create the SysManager object. This object will coordinate actions among various ports
	// and within nodes. In this example we use this object to setup and open our port.
	sFnd::SysManager* myMgr = sFnd::SysManager::Instance();                           //Create System Manager myMgr

	try
	{
		printf("SC-HUB SysManager is Initialized\n");   //Print to console
		sFnd::SysManager::FindComHubPorts(comHubPorts);
		printf("Found %d SC Hubs\n", comHubPorts.size());

		for (portCount = 0; portCount < comHubPorts.size() && portCount < NET_CONTROLLER_MAX; portCount++) {

			myMgr->ComHubPort(portCount, comHubPorts[portCount].c_str());    //define the first SC Hub port (port 0) to be associated 
											// with COM portnum (as seen in device manager)
		}

		if (portCount > 0) {
			//printf("\n I will now open port \t%i \n \n", portnum);
			myMgr->PortsOpen(portCount);             //Open the port


			sFnd::IPort& myPort = myMgr->Ports(0);

			printf(" Port[%d]: state=%d, nodes=%d\n",
				myPort.NetNumber(), myPort.OpenState(), myPort.NodeCount());

			sFnd::INode& theNodeX = myPort.Nodes(0);
			sFnd::INode& theNodeY = myPort.Nodes(1);
			sFnd::INode& theNodeZ = myPort.Nodes(2);

			theNodeX.EnableReq(false);
			theNodeY.EnableReq(false);
			theNodeZ.EnableReq(false);

			myMgr->Delay(200);
			theNodeX.Status.AlertsClear();
			theNodeX.Motion.NodeStopClear();
			theNodeX.EnableReq(true);

			printf("Node \t%zi enabled, Serial #: %d\n", 0, theNodeX.Info.SerialNumber.Value());
			double timeoutX = myMgr->TimeStampMsec() + TIME_TILL_TIMEOUT;

			theNodeY.Status.AlertsClear();
			theNodeY.Motion.NodeStopClear();
			theNodeY.EnableReq(true);

			printf("Node \t%zi enabled, Serial #: %d\n", 0, theNodeY.Info.SerialNumber.Value());
			double timeoutY = myMgr->TimeStampMsec() + TIME_TILL_TIMEOUT;

			theNodeZ.Status.AlertsClear();
			theNodeZ.Motion.NodeStopClear();
			theNodeZ.EnableReq(true);

			printf("Node \t%zi enabled, Serial #: %d\n", 0, theNodeZ.Info.SerialNumber.Value());
			double timeoutZ = myMgr->TimeStampMsec() + TIME_TILL_TIMEOUT;

			while (!theNodeX.Motion.IsReady()) {
				if (myMgr->TimeStampMsec() > timeoutX) {
					printf("Error: Timed out waiting for Node %d to enable\n", 0);
					msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
					return -2;
				}
			}

			while (!theNodeY.Motion.IsReady()) {
				if (myMgr->TimeStampMsec() > timeoutY) {
					printf("Error: Timed out waiting for Node %d to enable\n", 1);
					msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
					return -2;
				}
			}

			while (!theNodeZ.Motion.IsReady()) {
				if (myMgr->TimeStampMsec() > timeoutZ) {
					printf("Error: Timed out waiting for Node %d to enable\n", 2);
					msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
					return -2;
				}
			}

			theNodeX.Motion.MoveWentDone();
			theNodeX.AccUnit(sFnd::INode::RPM_PER_SEC);			    //Set the units for Acceleration to RPM/SEC
			theNodeX.VelUnit(sFnd::INode::RPM);					    //Set the units for Velocity to RPM
			theNodeX.Motion.AccLimit = ACC_LIM_RPM_PER_SEC *1;		//Set Acceleration Limit (RPM/Sec)
			theNodeX.Motion.VelLimit = VEL_LIM_RPM * 1; 		    //Set Velocity Limit (RPM)

			theNodeY.Motion.MoveWentDone();
			theNodeY.AccUnit(sFnd::INode::RPM_PER_SEC);			    //Set the units for Acceleration to RPM/SEC
			theNodeY.VelUnit(sFnd::INode::RPM);					    //Set the units for Velocity to RPM
			theNodeY.Motion.AccLimit = ACC_LIM_RPM_PER_SEC * 1;		//Set Acceleration Limit (RPM/Sec)
			theNodeY.Motion.VelLimit = VEL_LIM_RPM * 1; 		    //Set Velocity Limit (RPM)

			theNodeZ.Motion.MoveWentDone();
			theNodeZ.AccUnit(sFnd::INode::RPM_PER_SEC);			    //Set the units for Acceleration to RPM/SEC
			theNodeZ.VelUnit(sFnd::INode::RPM);					    //Set the units for Velocity to RPM
			theNodeZ.Motion.AccLimit = ACC_LIM_RPM_PER_SEC * 2;		//Set Acceleration Limit (RPM/Sec)
			theNodeZ.Motion.VelLimit = VEL_LIM_RPM * 2; 		    //Set Velocity Limit (RPM)

			size_t indexPoint = 1;
			for (const Vec3d& p : list)
			{
				Vec3d measurement;
				Vec3d background;

				Vec3d displacement = p - probePosition;
				
				int MOVE_DISTANCE_CNTSX = -std::round(displacement.x * xyCmToCounts);
				int MOVE_DISTANCE_CNTSY = std::round(displacement.y * xyCmToCounts);
				int MOVE_DISTANCE_CNTSZ = std::round(displacement.z * zCmToCounts);


				//printf("Moving Node \t%zi \n", 0);
				theNodeZ.Motion.MovePosnStart(MOVE_DISTANCE_CNTSZ);
				//printf("%f estimated time.\n", theNodeZ.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSZ));
				timeoutZ = myMgr->TimeStampMsec() + theNodeZ.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSZ) + 100;			//define a timeout in case the node is unable to enable

				while (!theNodeZ.Motion.MoveIsDone()) {
					if (myMgr->TimeStampMsec() > timeoutZ) {
						printf("Error: Timed out waiting for move to complete\n");
						msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
						return -2;
					}
				}
				//printf("Node \t%zi Move Done\n", 0);


				//printf("Moving Node \t%zi \n", 1);
				theNodeX.Motion.MovePosnStart(MOVE_DISTANCE_CNTSX);			//Execute 10000 encoder count move
				//printf("%f estimated time.\n", theNodeX.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSX));
				timeoutX = myMgr->TimeStampMsec() + theNodeX.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSX) + 100;			//define a timeout in case the node is unable to enable

				while (!theNodeX.Motion.MoveIsDone()) {
					if (myMgr->TimeStampMsec() > timeoutX) {
						printf("Error: Timed out waiting for move to complete\n");
						msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
						return -2;
					}
				}
				//printf("Node \t%zi Move Done\n", 1);


				//printf("Moving Node \t%zi \n", 2);
				theNodeY.Motion.MovePosnStart(MOVE_DISTANCE_CNTSY);			//Execute 10000 encoder count move
				//printf("%f estimated time.\n", theNodeY.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSY));
				timeoutY = myMgr->TimeStampMsec() + theNodeY.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSY) + 100;			//define a timeout in case the node is unable to enable

				while (!theNodeY.Motion.MoveIsDone()) {
					if (myMgr->TimeStampMsec() > timeoutY) {
						printf("Error: Timed out waiting for move to complete\n");
						msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
						return -2;
					}
				}
				//printf("Node \t%zi Move Done\n", 2);
				int sleepTime = 500;

				myMgr->Delay(1000);

				Sleep(5 * 1e3);

				if (back_flag == 0) {
					measurement = measureField(taskHandle, NumSamples, &read, data);
					std::cout << measurement << std::endl;
					float64 Bx = 0, By = 0, Bz = 0, V = 0;
					for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
						//write [t,Bx,By,Bz] to the file for each line
						Bx += data[i] / NumSamples;
						By += data[i + 1] / NumSamples;
						Bz += data[i + 2] / NumSamples;
						V += data[i + 3] / NumSamples;
					}
					float64 Bxsig = 0, Bysig = 0, Bzsig = 0, Vsig = 0;
					for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
						//write [t,Bx,By,Bz] to the file for each line
						Bxsig += (data[i] - Bx)* (data[i] - Bx) / NumSamples;
						Bysig += (data[i + 1]-By)* (data[i + 1] - By) / NumSamples;
						Bzsig += (data[i + 2]-Bz)* (data[i + 2] - Bz) / NumSamples;
						Vsig += (data[i + 3]-V)* (data[i + 3] - V) / NumSamples;
					}

					std::vector<double> outVals = { V, Bx, By, Bz, sqrt(Vsig), sqrt(Bxsig), sqrt(Bysig), sqrt(Bzsig)};

					std::cout << "Point #: " << indexPoint << "/" << list.size() << std::endl;
					indexPoint++;
					std::cout << "Bx : " << measurement.x << std::endl;

					//myMgr->Delay(500);
					//measurement.y = mm.measureBy();
					std::cout << "By : " << measurement.y << std::endl;

					//measurement.z = mm.measureBz();
					std::cout << "Bz : " << measurement.z << std::endl;
					//myMgr->Delay(500);

					result.push_back(outVals);
				}

				else if (back_flag == 1) {
					background = measureField(taskHandle, NumSamples, &read, data);

					float64 Bxb = 0, Byb = 0, Bzb = 0, Vb = 0;
					for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
						//write [t,Bx,By,Bz] to the file for each line
						Bxb += data[i] / NumSamples;
						Byb += data[i + 1] / NumSamples;
						Bzb += data[i + 2] / NumSamples;
						Vb += data[i + 3] / NumSamples;
					}
					float64 Bxbsig = 0, Bybsig = 0, Bzbsig = 0, Vbsig = 0;
					for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
						//write [t,Bx,By,Bz] to the file for each line
						Bxbsig += (data[i] - Bxb) * (data[i] - Bxb) / NumSamples;
						Bybsig += (data[i + 1] - Byb) * (data[i + 1] - Byb) / NumSamples;
						Bzbsig += (data[i + 2] - Bzb) * (data[i + 2] - Bzb) / NumSamples;
						Vbsig += (data[i + 3] - Vb) * (data[i + 3] - Vb) / NumSamples;
					}

					Status |= powerOn(vipwr, Status, TRUE);
					Sleep(5 * 1e3); 

					measurement = measureField(taskHandle, NumSamples, &read, data);
					float64 Bx = 0, By = 0, Bz = 0, V = 0;
					for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
						//write [t,Bx,By,Bz] to the file for each line
						Bx += data[i] / NumSamples;
						By += data[i + 1] / NumSamples;
						Bz += data[i + 2] / NumSamples;
						V += data[i + 3] / NumSamples;
					}
					float64 Bxsig = 0, Bysig = 0, Bzsig = 0, Vsig = 0;
					for (uInt32 i = 0; i < 4 * NumSamples; i += 4) {
						//write [t,Bx,By,Bz] to the file for each line
						Bxsig += (data[i] - Bx) * (data[i] - Bx) / NumSamples;
						Bysig += (data[i + 1] - By) * (data[i + 1] - By) / NumSamples;
						Bzsig += (data[i + 2] - Bz) * (data[i + 2] - Bz) / NumSamples;
						Vsig += (data[i + 3] - V) * (data[i + 3] - V) / NumSamples;
					}

					std::vector<double> outVals = { Vb, Bxb, Byb, Bzb, sqrt(Vbsig), sqrt(Bxbsig), sqrt(Bybsig), sqrt(Bzbsig), V, Bx, By, Bz, sqrt(Vsig), sqrt(Bxsig), sqrt(Bysig), sqrt(Bzsig)};
					Status |= powerOn(vipwr, Status, FALSE);

					result.push_back(outVals);

				}
				

				probePosition = probePosition + displacement;
			}

			Vec3d origin{ 0,0,0 };
			Vec3d displacement = origin - probePosition;

			int MOVE_DISTANCE_CNTSZ = std::round(displacement.z * zCmToCounts);
			int MOVE_DISTANCE_CNTSX = -std::round(displacement.x * xyCmToCounts);
			int MOVE_DISTANCE_CNTSY = std::round(displacement.y * xyCmToCounts);

			//printf("Moving Node \t%zi \n", 0);
			theNodeZ.Motion.MovePosnStart(MOVE_DISTANCE_CNTSZ);			//Execute 10000 encoder count move
			//printf("%f estimated time.\n", theNodeZ.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSZ));
			timeoutZ = myMgr->TimeStampMsec() + theNodeZ.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSZ) + 100;			//define a timeout in case the node is unable to enable

			while (!theNodeZ.Motion.MoveIsDone()) {
				if (myMgr->TimeStampMsec() > timeoutZ) {
					printf("Error: Timed out waiting for move to complete\n");
					msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
					return -2;
				}
			}
			//printf("Node \t%zi Move Done\n", 0);


			//printf("Moving Node \t%zi \n", 1);
			theNodeX.Motion.MovePosnStart(MOVE_DISTANCE_CNTSX);			//Execute 10000 encoder count move
			//printf("%f estimated time.\n", theNodeX.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSX));
			timeoutX = myMgr->TimeStampMsec() + theNodeX.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSX) + 100;			//define a timeout in case the node is unable to enable

			while (!theNodeX.Motion.MoveIsDone()) {
				if (myMgr->TimeStampMsec() > timeoutX) {
					printf("Error: Timed out waiting for move to complete\n");
					msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
					return -2;
				}
			}
			//printf("Node \t%zi Move Done\n", 1);


			//printf("Moving Node \t%zi \n", 2);
			theNodeY.Motion.MovePosnStart(MOVE_DISTANCE_CNTSY);			//Execute 10000 encoder count move
			//printf("%f estimated time.\n", theNodeY.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSY));
			timeoutY = myMgr->TimeStampMsec() + theNodeY.Motion.MovePosnDurationMsec(MOVE_DISTANCE_CNTSY) + 100;			//define a timeout in case the node is unable to enable

			while (!theNodeY.Motion.MoveIsDone()) {
				if (myMgr->TimeStampMsec() > timeoutY) {
					printf("Error: Timed out waiting for move to complete\n");
					msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key
					return -2;
				}
			}
			//printf("Node \t%zi Move Done\n", 2);

			// Disable all nodes
			myPort.Nodes(0).EnableReq(false);
			myPort.Nodes(1).EnableReq(false);
			myPort.Nodes(2).EnableReq(false);

		}
		else {
			printf("Unable to locate SC hub port\n");

			msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key

			return -1;  //This terminates the main program
		}

	}
	catch (sFnd::mnErr & theErr)    //This catch statement will intercept any error from the Class library
	{
		printf("Port Failed to open, Check to ensure correct Port number and that ClearView is not using the Port\n");
		//This statement will print the address of the error, the error code (defined by the mnErr class), 
		//as well as the corresponding error message.
		printf("Caught error: addr=%d, err=0x%08x\nmsg=%s\n", theErr.TheAddr, theErr.ErrorCode, theErr.ErrorMsg);

		msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key

		return -1;  //This terminates the main program
	}

	myMgr->PortsClose();
	delete[] data;

}

//ViStatus powerOn(ViSession& pwr, ViStatus Status, bool on) {
//	if (on) {
//		Status |= viPrintf(pwr, "OUTP ON\n");
//	}
//	else {
//		Status |= viPrintf(pwr, "OUTP OFF\n");
//	}
//	return Status; 
//}

//ViStatus powerOn(ViSession& pwr, ViStatus Status) {
//	Status |= viPrintf(pwr, "OUTP ON\n");
//	return Status;
//}

//void setCurrent(ViSession& pwr, ViStatus Status, float curr) {
//	Status |= viPrintf(pwr, "CURR %.4f\n", curr);
//}

int main()
{
	std::cout << " ---------------------------------------------" << std::endl;
	std::cout << "|  Command Line Auto Field Mapper v2.0.1      |" << std::endl;
	std::cout << "|  Umit H. Coskun. Last update : 2/20/2020    |" << std::endl;
	std::cout << "|  Andrew Mullins  Last update : 11/21/2023   |" << std::endl;
	std::cout << "|  Richard McDonald  Last update : 3/28/2025  |" << std::endl;
	std::cout << " ---------------------------------------------" << std::endl;

	std::ofstream outfile;
	// setup the NIDAQmx object
	int32           error = 0;
	TaskHandle      taskHandle = 0;
	int32           read;
	char            errBuff[2048] = { '\0' };
	uInt32          SampleRate = 5000; // in Hz
	uInt32          NumSamples = 5000;
	std::string     fileName = "dataout.txt";

	int32			back_flag = 0;

	MatrixXd w = MatrixXd::Zero(NumSamples,3);

	// Create the read task to read from the 3 channels (0-2)
	DAQmxCreateTask("", &taskHandle);
	DAQmxCreateAIVoltageChan(taskHandle, "cDAQ1Mod1/ai0:3", "", DAQmx_Val_Cfg_Default, -10.0, 10.0, DAQmx_Val_Volts, NULL);
	DAQmxCfgSampClkTiming(taskHandle, "", SampleRate, DAQmx_Val_Rising, DAQmx_Val_ContSamps, NumSamples);

	// DAQmx Start Code
	DAQmxStartTask(taskHandle);

	// Power Supply 
	const char pwrVISAAddress[14] = "ASRL14::INSTR";  // power supply ; might need to check NI MAX for port #
	ViConstRsrc pwr = pwrVISAAddress;
	ViStatus Status;
	ViSession defaultRMpwr, vipwr;
	Status = viOpenDefaultRM(&defaultRMpwr);
	Status |= viOpen(defaultRMpwr, pwr, VI_NULL, VI_NULL, &vipwr);

	// Declare variables
	std::vector<Vec3d> coordinateList;
	std::string outputPath;
	std::vector<std::vector<double>> resultList;
	char input;
	//dmm.getVISA_Addr();
	while (1)
	{
		std::cout << " ----------------------------------------------------------------" << std::endl;
		std::cout << "|  Test measure[t], Load coordinates[l], Print coordinates[p]     |" << std::endl;
		std::cout << "|  Enter output file name[o], Start Mapping[m], Power ON[b],      |" << std::endl;
		std::cout << "|  Power OFF[c], Set Current[i], Toggle Background Subtraction[g] |" << std::endl;
		std::cout << "|  Quit[q]														|" << std::endl;
		std::cout << " ---------------------------------------------------------------- " << std::endl;
		std::cout << "Command[]: ";

		std::cout << std::setw(1);
		std::cin >> input;
			switch (input)
			{
			case 'b':
			{
				Status |= powerOn(vipwr, Status, TRUE);
			}
			break;
			case 'c':
			{
				Status |= powerOn(vipwr, Status, FALSE);
			}
			break;

			case 'i':
			{
				float curr;
				std::cout << "Enter the current in Amps: ";
				std::cin >> curr;
				setCurrent(vipwr, Status, curr);
			}
			break;

			case 'g':
			{
				isSubtractionOn = !isSubtractionOn;
				std::cout << "Background Subtraction : " << isSubtractionOn << std::endl;
				back_flag = 1;
			}
			break;

			case 't':
			{
				testMeasureField(taskHandle, NumSamples);
			}
			break;

			case 'q':
			{
				return 0;
			}
			break;

			case 'p':
			{
				printVec3dList(coordinateList);
			}
			break;

			case 'l':
			{
				std::string path;
				std::cout << "Enter the coordinate file name: ";
				std::cin >> path;
				FileReader reader(path);
				coordinateList = reader.getCoordinates();
			}
			break;

			case 'o':
			{
				std::cout << "Enter the output file name: ";
				std::cin >> outputPath;
				std::cout << "Output filename set to : " << outputPath << std::endl;
			}
			break;

			case 'm':
			{
				std::cout << "This will start automatic mapping. Please verify following: " << std::endl;
				std::cout << "Coordinates :" << std::endl;
				printVec3dList(coordinateList);
				std::cout << "Output file name: " << std::endl;
				std::cout << outputPath << std::endl;

				//   Below is AC current handling written by Andrew 
				// no error handling here, so just don't put in something that isn't a double...
				std::cout << "What is the measurement frequency?: ";
				double mFreq;
				std::cin >> mFreq;

				// We are going to make the pinv matrix here since it will be the same for every point. 
				// We have to wait until the user provides the input frequency before we can generate the array.
				VectorXd t = VectorXd::LinSpaced(NumSamples, 0, NumSamples / SampleRate);
				w.col(0) = (2 * M_PI * t * mFreq).array().sin();
				w.col(1) = (2 * M_PI * t * mFreq).array().cos();
				w.col(2) = VectorXd::Ones(NumSamples);

				MatrixXd wt = w.transpose();
				MatrixXd pinv = (wt*w).inverse() * wt;

				std::cout << "Do you really want to start mapping ? [y/n]: ";
				char answer;
				std::cin >> answer;
				switch (answer)
				{
				case 'y':
				{
					if (back_flag == 0) {
						outfile.open(outputPath);
						Map(coordinateList, taskHandle, NumSamples, resultList, pinv, back_flag, vipwr, Status);
						//printVec3dList(resultList);
						for (size_t i = 0; i < coordinateList.size(); i++) {
							std::string outString = "";

							for (int j = 0; j < resultList[i].size(); j++) {
								outString = outString + " " + std::to_string(resultList[i][j]);
							}
							outfile << coordinateList.at(i).x << " " << coordinateList.at(i).y << " " << coordinateList.at(i).z << " " << outString << std::endl;

						}
						outfile.close();
					}
					else if (back_flag == 1) {
						outfile.open(outputPath);
						Map(coordinateList, taskHandle, NumSamples, resultList, pinv, back_flag, vipwr, Status);
						std::cout << coordinateList.size() << std::endl;
						std::cout << resultList.size() << std::endl;
						//printVec3dList(resultList);
						for (size_t i = 0; i < coordinateList.size(); i++) {
							std::string outString = "";

							for (int j = 0; j < resultList[i].size(); j++) {
								outString = outString + " " + std::to_string(resultList[i][j]);
							}
							outfile << coordinateList.at(i).x << " " << coordinateList.at(i).y << " " << coordinateList.at(i).z << outString << std::endl;

						}
						outfile.close();
					}
				}
				break;
				}
			}
			break;
			}
		
	}
}
