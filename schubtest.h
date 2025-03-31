#pragma once
#include "geometry.h"
#include <stdio.h>    
#include <iostream>
#include "pubSysCls.h"  

// Send message and wait for newline
void msgUser(const char* msg) {
	std::cout << msg;
	getchar();
}

class MotorControl
{
public:
	int xyCmToCounts = 750;
	int zCmToCounts = 3000;

	int ACC_LIM_RPM_PER_SEC = 20;
	int VEL_LIM_RPM = 50;
	int TIME_TILL_TIMEOUT = 10000;
	int MOVE_DISTANCE_CNTS;

	size_t portCount = 0;
	std::vector<std::string> comHubPorts;
	//Create the SysManager object. This object will coordinate actions among various ports
	// and within nodes. In this example we use this object to setup and open our port.
	sFnd::SysManager* myMgr = sFnd::SysManager::Instance();                           //Create System Manager myMgr
	MotorControl()
	{
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

				static sFnd::INode& theNodeZ = myPort.Nodes(0);
				theNodeZ.EnableReq(false);

				//static sFnd::INode& theNodeX = myPort.Nodes(1);
				//static sFnd::INode& theNodeY = myPort.Nodes(2);

				printf(" Port[%d]: state=%d, nodes=%d\n",
					myPort.NetNumber(), myPort.OpenState(), myPort.NodeCount());

			}
			else {
				printf("Unable to locate SC hub port\n");

				msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key

				return;  //This terminates the main program
			}

		}
		catch (sFnd::mnErr & theErr)    //This catch statement will intercept any error from the Class library
		{
			printf("Port Failed to open, Check to ensure correct Port number and that ClearView is not using the Port\n");
			//This statement will print the address of the error, the error code (defined by the mnErr class), 
			//as well as the corresponding error message.
			printf("Caught error: addr=%d, err=0x%08x\nmsg=%s\n", theErr.TheAddr, theErr.ErrorCode, theErr.ErrorMsg);

			msgUser("Press any key to continue."); //pause so the user can see the error message; waits for user to press a key

			return;  //This terminates the main program
		}

	}

	~MotorControl()
	{
		// Close down the ports
		myMgr->PortsClose();
	}

};
