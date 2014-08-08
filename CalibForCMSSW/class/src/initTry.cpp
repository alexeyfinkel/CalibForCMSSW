//Here, I try to implement some of the basic funcitons for the calibration class; still in brainstorming/trial stage.


////constructor////
EECalibration::EECalibration(TChain inputChain, std::string outputDir)
{
	initHists();
	initMaps();
	if(histsReady && mapsReady)
	{
		calibrate();
	}
	pruneHists();
	printConsts();
}
//====================================================================

////Once per Job////
EECalibration::calibrate()
{
	for(int iteration=1;iteration<=nIterations;iteration++)
	{
		initTree();//return to start of tree if iteration =/= 1
		if(treeReady)
		{
			fillRatHists();
		}
		updateConsts();
		//***remove tree***//
	}
}

void EECalibration::initHists()
{
	calConstsN = newTH2D("calConstsN","Calibration Constants, EE-;i_{x};i_{y}",100,0,100);
	calConstsP = newTH2D("calConstsP","Calibration Constants, EE+;i_{x};i_{y}",100,0,100);
	
	for(int i=1;i<101;i++)//doing all of EE now! (Woo!)
	{
		for(int j=1;j<101;j++)
		{
			calConstsN->SetBinContent(i,j,1);
			calConstsP->SetBinContent(i,j,1);
			calConstsN->SetBinError(i,j,999);
			calConstsP->SetBinError(i,j,999);
		}
	}
	histsReady = true;
}

void EECalibration::initMaps()
{
	char name[128];
	for(int i=1;i<101;i++)//doing all of EE now! (Woo!)
	{
		for(int j=1;j<101;j++)
		{
			sprintf(name,"ratioN_%d_%d",i,j);
			ratioMapN[std::make_pair(i,j)] = new TH1D(name,"",24,0.4,1.6);
			ratioMapN[std::make_pair(i,j)]->SetDefaultSumw2();//since will be using weights =/=1
			sprintf(name,"ratioP_%d_%d",i,j);
			ratioMapP[std::make_pair(i,j)] = new TH1D(name,"",24,0.4,1.6);
			ratioMapP[std::make_pair(i,j)]->SetDefaultSumw2();
		}
	}
	mapsReady = true;
}

void EECalibration::pruneHists()
{
	for(int i=1;i<101;i++)
	{
		for(int j=1;j<101;j++)
		{
			if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 2550.25) )
			{
				if(constsN->GetBinContent(i,j)<0.4 || constsN->GetBinContent(i,j)>1.6 || constsN->GetBinError(i,j)>5)
				{
					constsN->SetBinContent(i,j,-1);
					constsN->SetBinError(i,j,999);
				}
				if(constsP->GetBinContent(i,j)<0.4 || constsN->GetBinContent(i,j)>1.6 || constsP->GetBinError(i,j)>5)
				{
					constsP->SetBinContent(i,j,-1);
					constsP->SetBinError(i,j,999);
				}
			}
		}
	}
}

void EECalibration::printConsts()
{
	char fname[256];
	sprintf(fname,"%s/ZeeCalibrationConstants.txt",outputDir.c_str());
	FILE* constsFile;
	constsFile = fopen(fname,"w");
	for(int i=1;i<101;i++)
	{
		for(int j=1;j<101;j++)
		{
			if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 2550.25) )
			{
				fprintf(constsFile,"%d\t%d\t%d\t%.6f\t%.6f\n",i,j,-1,constsN->GetBinContent(i,j),GetBinError(i,j));
				fprintf(constsFile,"%d\t%d\t%d\t%.6f\t%.6f\n",i,j,-1,constsN->GetBinContent(i,j),GetBinError(i,j));
			}
		}
	}
}
//==============================================================================================================================

////Once Per Iteration////
void EECalibration::initTree()
{
	//turn off all branches, and reactivate only the useful ones:
	inputChain->SetBranchStatus("*",0);
	inputChain->SetBranchStatus("XRecHitSCEle",1);//activate hix
	inputChain->SetBranchStatus("YRecHitSCEle",1);//activate hiy
	inputChain->SetBranchStatus("",1);//activate fractions
	inputChain->SetBranchStatus("",1);//activate expected energy
	inputChain->SetBranchStatus(EName.c_str(),1);//activate observed energy of your choosing
	
	//set branch addresses:
	inputChain->SetBranchAddress("", &expectedE);//***expectedE branch forthcoming...
	inputChain->SetBranchAddress(EName.c_str(), observedEs);
	inputChain->SetBranchAddress("XRecHitSCEle", hix);
	inputChain->SetBranchAddress("YRecHitSCEle", hiy);
	inputChain->SetBranchAddress("", hitFractions); //***fractions branch forthcoming...
	
	inputChain->GetEvent(entryList->GetEntry(0));
}

void EECalibration::resetRatHists()
{
	for(int i=1;i<101;i++)
	{
		for(int j=1;j<101;j++)
		{
			ratioMapP[std::male_pair(i,j)]->Reset();
			ratioMapP[std::male_pair(i,j)]->Reset();
			
		}
	}
}

void EECalibration::fillRatHists()
{
	for(int evt=0;evt<entryList->GetN();evt++)
	{
		clearVars(i);//***This needs to be updated... once the final vars content is decided***//
		if(!getElectrons(entryList->GetEntry(evt)))) continue;//try to get electrons in correct range, passing selection
		if(iteration==1) correction = 1;
		else
		{
			for(unsigned int hit=0; hit<hitfractions[probe].size_of();hit++)
			{
				if(etas[probe]<0)
					correction += hitFractions[probe].at(hit)*constsN->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) );
					
				else correction += hitFractions[probe].at(hit)*constsP->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) );
			}
		}
		energyRatio = expectedE/(observedEs[probe]*correction);
		
		for(unsigned int hit=0; hit<hitEnergies[probe].size_of();hit++)
		{
			if(etas[probe]<0)
				ratioMapN[std::male_pair( hix[probe].at(hit),hiy[probe].at(hit) )]->Fill( energyRatio,hitFractions[probe].at(hit)
							*constsN->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) );
				
			else ratioMapP[std::male_pair( hix[probe].at(hit),hiy[probe].at(hit) )]->Fill( energyRatio,hitFractions[probe].at(hit)
							*constsP->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) );
		}		
	}	
}

void EECalibration::updateConsts()
{
	for(int i=1;i<101;i++)
	{
		for(int j=1;j<101;j++)
		{
			if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 2550.25) )
			{
				//ze2->GetNextEvent(); need correct version...
				continue;
			}
			xtlFit(ratioMapN[std::male_pair(i,j)]);
			constsN->SetBinContent(i,j, constsN->GetBinContent(i,j)*gaussRatioFit->GetParameter(1) );
			constsN->SetBinError(i,j, gaussRatioFit->GetParError(1) );
			
			xtlFit(ratioMapP[std::male_pair(i,j)]);//internal safety checks should catch bad fits now
			constsP->SetBinContent(i,j, constsN->GetBinContent(i,j)*gaussRatioFit->GetParameter(1) );
			constsP->SetBinError(i,j, gaussRatioFit->GetParError(1) );
		}
	}
}

void EECalibration::clearVars()//***Will need to be updated once needed variables are settled down!...***//
{
	expectedE = 0;
	observedEs.clear();
	hitFractions.clear();
	hix.clear(); //0 should be the safest value, as it should correspond to underflow bins in cal. const. hists
	hiy.clear();
}
//===============================================================================================================================

////Once Per Event, or Many Times Per Iteration
bool EECalibration::getElectrons(int i)
{
	inputChain->GetEvent(i);
	//sort electrons by location
	if(etas[0]<1.48 && etas[1]>1.57)
	{
		tag = 0;
		probe = 1;
	}
	else if(etas[1]<1.48 && etas[0]>1.57)
	{
		tag = 1;
		probe = 0;
	}
	else return false;
	
	//assuming range is fine:	
	return true;		
}

void EECalibration::xtlFit(TH1D* xtlHist)
{
	//prepare reasonable parameter values
	gfit->SetParameter(0,10);
	gfit->SetParameter(1,1);
	gfit->SetParLimits(1,0.2,1.8);
	gfit->SetParameter(2,0.1);
	gfit->SetParLimits(2,0,1);
	//fit...
	xtlHist->Fit(gfit,"QLMN",0.4,1.6);
	//check for runaways:
	if(gfit->GetParameter(1)==0.2 || gfit->GetParameter(1)==1.8) gfit->SetParameter(1,-1);
	if(gfit->GetParameter(2)==0 || gfit->GetParameter(2)==1) gfit->SetParameter(1,999);
}















