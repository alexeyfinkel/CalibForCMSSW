//Here, I try to implement some of the basic funcitons for the calibration class; still in brainstorming/trial stage.



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

void EECalibration::initTree()
{
	
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
	for(int evt=0;evt<inputChain->GetEntries();evt++)
	{
		initVars();//***This needs to be updated... once the final vars content is decided***//
		if(!getElectrons(i)) continue;//try to get electrons in correct range, passing selection
		expectedPt = Mz*Mz/( pts[tag]*cosh(etas[tag]-etas[probe])-cos(phis[tag]-phis[probe]) );
		if(iteration==1) correction = 1;
		else
		{
			for(unsigned int hit=0; hit<hitEnergies[probe].size_of();hit++)
			{
				if(etas[probe]<0)
				{
					correctedE += hitEnergies[probe].at(hit)*constsN->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) );
					hitFractions.push_back( hitEnergies[probe].at(hit)*constsN->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) )/rawEnergies[probe].at(hit) );
				}
				else 
				{	
					correctedE += hitEnergies[probe].at(hit)*constsP->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) );
					hitFractions.push_back( hitEnergies[probe].at(hit)*constsP->GetBinContent( hix[probe].at(hit),hiy[probe].at(hit) )/rawEnergies[probe].at(hit) );
				
			}
			correction = correctedE/rawE;
		}
		energyRatio = expectedPt/(pts[probe]*correction);
		
		for(unsigned int hit=0; hit<hitEnergies[probe].size_of();hit++)
		{
			if(etas[probe]<0)
				ratioMapN[std::male_pair( hix[probe].at(hit),hiy[probe].at(hit) )]->Fill( energyRatio,hitFractionsat(hit) );
			else ratioMapP[std::male_pair( hix[probe].at(hit),hiy[probe].at(hit) )]->Fill( energyRatio,hitFractionsat(hit) );
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

void EECalibration::initVars()//***Will need to be updated once needed variables are settled down!...***//
{
	Z.e1_pt = -1000;
	Z.e1_eta = -100;
	Z.e1_phi = -100;
	Z.e2_pt = -1000;
	Z.e2_eta = -100;
	Z.e2_phi = -100;
	Z.hits.clear();
	Z.eFractions.clear();
	hix = 0; //0 should be the safest value, as it should correspond to underflow bins in cal. const. hists
	hiy = 0;
}

bool EECalibration::getElectrons(int evt)
{
	inputChain->GetEvent(evt);
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
	//***Need to either check appropriate selection bit or manually check selection criteria
	//This is different for Tracked and NT electrons!***//
	//assuming electrons pass appropriate selection:
	
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















