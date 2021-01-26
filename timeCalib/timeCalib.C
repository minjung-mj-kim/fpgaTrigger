#include <iostream>
#include <string>
#include <cstdlib>
#include <tuple>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TString.h>
#include <TObjArray.h>

/** Trig numbers **/
#define N_TRIG 12  
#define N_BOARD 5 
#define N_CH 96 
#define N_DISC 6 // (5 discriminators +  no discriminator)
#define N_STATION 4

/** TDC unit to nano second conversion **/
#define TDC_NS_CONV 1// 1 is for no conversion, 25/16 for 40 MHz clock

/** Reference time for delay adjustment **/
#define TIME_REF 650

/** Hit multiplicity stuff (to check noise level in the future diagnosis, not required for timing calib) **/
#define MAX_MULTI 10 

/** TDC range to analyze **/
#define TDC_MIN 500
#define TDC_MAX 650

/** Plot stuffs **/
const float TIME_MIN = (TDC_MIN*TDC_NS_CONV);
const float TIME_MAX = (TDC_MAX*TDC_NS_CONV);
const int   TIME_BIN = (TDC_MAX-TDC_MIN);


using namespace std;


typedef tuple<int, string, string> tuple_info;
typedef tuple<int, string, string, string, string> tuple_v1495;

const tuple_info info_trig[N_TRIG] =
  {
   tuple_info(0, "FPGA1", "ST 2 & 4"),
   tuple_info(1, "FPGA2", "ST 2T & 4T"),
   tuple_info(2, "FPGA3", "ST 2T & 4B"),
   tuple_info(3, "FPGA4", "ST 2B & 4T"),
   //tuple_info(4, "FPGA5", "NIM ST 2 & 4"), // ~run 1501
   tuple_info(4, "FPGA5", "ST 2B & 4B"), // run 1502~  
   tuple_info(5, "NIM1", "ST 1 & 2 & 3 & 4"),
   tuple_info(6, "NIM2", "ST 1 & 2"),
   tuple_info(7, "NIM3", "Random"),
   // tuple_info(8, "NIM4", "ST 2 & 3"), // ~run 1501
   tuple_info(8, "NIM4", "ST 2 & 4"), // run 1502~   
   tuple_info(9, "NIM5", "Flush"),
   tuple_info(10, "BOS", "Begin of spill"),
   tuple_info(11, "EOS", "End of spill")   
  };

/** Hodoscope Discriminator **/
const tuple_info info_disc[N_DISC] =
  {
   tuple_info(0, "ST1" , "ST1 H5-19"),
   tuple_info(1, "ST2" , "ST2 H1-16"),
   tuple_info(2, "ST3" , "ST3 H1-16"),
   tuple_info(3, "ST4a", "ST4 H1-8" ),
   tuple_info(4, "ST4b", "ST4 H9-16"),
   tuple_info(5, "NotHodo", "G port or Empty channel")
  };


/** delay and mapping files

    up-to-date timing files are under
    /home/e1039daq/noahKnowsBest/vme_workdir

    mapping files from redmine
**/
const tuple_v1495 info_board[N_BOARD] =
  {
   tuple_v1495(0, "420", "XT Lv-A", "Mapping/420mapping.txt", "../Timing/time_0.txt"),
   tuple_v1495(1, "430", "XB Lv-A", "Mapping/430mapping.txt", "../Timing/time_1.txt"),
   tuple_v1495(2, "460", "XT Lv-B", "Mapping/460mapping.txt", "../Timing/time_2.txt"),
   tuple_v1495(3, "470", "XB Lv-B", "Mapping/470mapping.txt", "../Timing/time_3.txt"),
   tuple_v1495(4, "480", "XT Lv-C", "Mapping/480mapping.txt", "../Timing/time_4.txt")
  };

struct TDC
{
  string ch; // v1495 TDC channel
  string disc; // discriminator module
  string board; // v1495 board
  string station; // hodoscope station
  string port; // v1495 port name (A,B,C,D,E,F,G), half port share a ribbon cable from discriminator
  string vhdl; // vhdl name. not used for now
  string name; // whole name for cross check
};


/** fit histogram to get timing peak 
    currently get peak bean for non gausian cosmic timings
    it will be changed to fit gaus for beam
**/
float get_peak(TH1* h)
{
  int peak = h->GetBinCenter(h->GetMaximumBin());
  return peak;
}

/** to make new delay file based on timing fit **/
int update_delay(float old_delay, float peak)
{
  int tdc_peak = int(peak/TDC_NS_CONV);
  int new_delay = (tdc_peak-TIME_REF)+old_delay;

  return new_delay;
}


int main(int argc, char* argv[])
{
  if( argc<3 )
    {
      cout << "Usage: ./timeCalib "
           << "arg1 = <output file tag, like runnumber or anything you want to name> " << endl      
           << "arg2, arg3, ... = <input runnumber1> <runnumber2>  ..." << endl;
      return 0;
    }

  
  // ====================================================
  //
  // pass argument (get runnumber, set file name)
  //
  // ====================================================

  
  //TFile* wfile = new TFile(argv[1],"RECREATE");
  
  const int N_FILE = (argc-2);
  TString rfile_name[N_FILE];
  TString run_name = "run";
  
  for(int irf=0; irf<N_FILE; irf++)
    {
      TString runnumber(argv[irf+2]);
      rfile_name[irf] = "datafile/run_" + runnumber + ".root";
      run_name += " " + runnumber;
    }

  
  // ====================================================
  //
  // Fill TDC info
  //
  // ====================================================

  
  TDC tdc[N_BOARD][N_CH];

  for(int iboard=0; iboard<N_BOARD; iboard++)
    {
      ifstream ifs(get<3>(info_board[iboard]));            

      TString line_chmap;
      line_chmap.ReadLine(ifs);

      
      for(int ich=0; ich<N_CH; ich++)
        {
          line_chmap.ReadLine(ifs);

          TObjArray * tempArray = line_chmap.Tokenize(",");

          TDC* ptr_tdc = &tdc[iboard][ich];

          ptr_tdc->ch    = tempArray->At(0)->GetName();
          ptr_tdc->name  = tempArray->At(1)->GetName();
          ptr_tdc->port  = tempArray->At(2)->GetName();
          ptr_tdc->vhdl  = tempArray->At(3)->GetName();
          ptr_tdc->board = get<1>(info_board[iboard]);
          
         
          if( ich < 16 )
            {
              ptr_tdc->disc = "3";
              ptr_tdc->station = "4";
            }
          else if( ich < 32 )
            {
              ptr_tdc->disc = "4";
              ptr_tdc->station = "4";
            }
          else if( ich < 48 )
            {
              ptr_tdc->disc = "1";
              ptr_tdc->station = "2";
            }
          else if( ich < 64 )
            {
              ptr_tdc->disc = "2";
              ptr_tdc->station = "3";
            }
          else if( ich > 65 && ich <80 )
            {              
              ptr_tdc->disc = "0";
              ptr_tdc->station = "1";
            }
          else if( ich == 64 )
            {              
              ptr_tdc->disc = "5";
              ptr_tdc->station = "5";
            }
          else if( ich == 65 )
            {              
              ptr_tdc->disc = "5";
              ptr_tdc->station = "6";
            }
          else
            {              
              ptr_tdc->disc = "5";
              ptr_tdc->station = "7";
            }

          if( ich!= stoi(ptr_tdc->ch) )
            cout << "scary bug" << endl;
        }
    }

  cout << "TDC mapping check: board, channel, name, station/else, disc, v1495port, vhdl " << endl;
  
  for(int iboard=0; iboard<N_BOARD; iboard++)
    {
      for(int ich=0; ich<N_CH; ich++)
        {
          cout << tdc[iboard][ich].board << "\t"
               << tdc[iboard][ich].ch << "\t"
               << tdc[iboard][ich].name << "\t"
               << tdc[iboard][ich].station << "\t"
               << tdc[iboard][ich].disc << "\t"
               << tdc[iboard][ich].port << "\t"
               << tdc[iboard][ich].vhdl << endl;          
        }
     }

  
  // ====================================================
  //
  // Initialize Histogram 
  //
  // "hh_tdc" is for presentation of timing calibration.
  //
  // The rest are for diagnosis, for paranoid people like me.
  //
  // ====================================================
  
  TH2F* hh_tdc   [N_TRIG][N_BOARD]; // time vs chID
  
  TH1F* h_tdc    [N_TRIG][N_BOARD][N_CH]; // time
  TH2F* hh_t_m   [N_TRIG][N_BOARD][N_CH]; // time vs #multiplicity

  TH1F* h_hit    [N_TRIG][N_BOARD]; // #hit vs chID
  TH2F* hh_multi [N_TRIG][N_BOARD]; // multiplicity vs chID
  
  TH1F* h_multi     [N_TRIG][N_BOARD]; // #multiplicity per board 
  TH1F* h_multi_st  [N_TRIG][N_BOARD][N_STATION]; // #multiplicity per station per board
  TH1F* h_multi_hst [N_TRIG][N_STATION]; // #hodo multiplicity per station
  TH1F* h_multi_all [N_TRIG]; // #hodo multiplicity per event
  
  
  
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      string obj_id;
      string label;

      obj_id = "hMultiAll_trig" + to_string(itrig);
      label  = run_name + " " + get<1>(info_trig[itrig]) + " ; Hodo hit multiplicity;  #event";      
      
      h_multi_all[itrig] = new TH1F( obj_id.c_str(), label.c_str(),
                                     MAX_MULTI*200, 0, MAX_MULTI*200);
      
      for(int istation=0; istation<N_STATION; istation++)
        {
          obj_id = "hMultiHodoSt_trig"+to_string(itrig)+"_st"+to_string(istation+1);
          label = run_name + " " + get<1>(info_trig[itrig]) + " ST"+to_string(istation+1)+ " ; #TDC hit per event; #event";
          
          h_multi_hst[itrig][istation] = new TH1F( obj_id.c_str(), label.c_str(),
                                                          MAX_MULTI*40, 0, MAX_MULTI*40);
        }
   
      for(int iboard=0; iboard<N_BOARD; iboard++)
        {
          obj_id = "hHit_trig" + to_string(itrig) + "_board" + to_string(iboard);
          label  = run_name + " " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard])+ "; Channel ID;  hit count";      
          
          h_hit[itrig][iboard] = new TH1F( obj_id.c_str(), label.c_str(),
                                           N_CH, 0, N_CH);
          
          
          obj_id = "hhTdc_trig" + to_string(itrig) + "_board" + to_string(iboard);
          label  = run_name + " " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard])+ "; Channel ID; Time (ns)";      
          
          hh_tdc[itrig][iboard] = new TH2F( obj_id.c_str(), label.c_str(),
                                            N_CH, 0, N_CH,
                                            TIME_BIN, TIME_MIN, TIME_MAX);
          
          
          obj_id = "hhMulti_trig" + to_string(itrig) + "_board" + to_string(iboard);
          label  = run_name + " " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard])+ "; Channel ID ; Multiplicity";      
          
          hh_multi[itrig][iboard] = new TH2F( obj_id.c_str(), label.c_str(),
                                              N_CH, 0, N_CH,
                                              MAX_MULTI, 0, MAX_MULTI);

          
          obj_id = "hMulti_trig"+to_string(itrig)+"_board"+to_string(iboard);
          label = run_name + " " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard]) + " ; #TDC hit per event; #event";
          
          h_multi[itrig][iboard] = new TH1F( obj_id.c_str(), label.c_str(),
                                             MAX_MULTI*100, 0, MAX_MULTI*100);


          for(int istation=0; istation<N_STATION; istation++)
            {
              obj_id = "hMultiSt_trig"+to_string(itrig)+"_board"+to_string(iboard)+"_st"+to_string(istation+1);
              label = run_name + " " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard]) + " ST"+to_string(istation+1)+ " ; #TDC hit per event; #event";
          
              h_multi_st[itrig][iboard][istation] = new TH1F( obj_id.c_str(), label.c_str(),
                                             MAX_MULTI*20, 0, MAX_MULTI*20);
            }
            
          for(int ich=0; ich<N_CH; ich++)
            {            
              obj_id = "hTdc_trig"+to_string(itrig)+"_board"+to_string(iboard)+"_ch"+to_string(ich);
              label = run_name + " " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard]) + " ch"+ich + "; Time (ns) ; #hit";
            
              h_tdc[itrig][iboard][ich] = new TH1F( obj_id.c_str(), label.c_str(),
                                                    TIME_BIN, TIME_MIN, TIME_MAX);
              
              obj_id = "hhTM_trig"+to_string(itrig)+"_board"+to_string(iboard)+"_ch"+to_string(ich);
              label = run_name + " " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard]) + " ch"+ich + " ; Multiplicity; Time (ns)";
              
              hh_t_m[itrig][iboard][ich] = new TH2F( obj_id.c_str(), label.c_str(),
                                                     MAX_MULTI, 0, MAX_MULTI,
                                                     TIME_BIN, TIME_MIN, TIME_MAX);
            }
        }
    }
  
  
  // ====================================================
  //
  // Event loop -> Fill histogram
  //
  // ====================================================

  
  unsigned eventID;
  unsigned eventType;
  unsigned triggerType;
  unsigned nHits;
  unsigned boardID[1000];
  unsigned channelID[1000];
  unsigned tdcTime[1000];
  unsigned triggerTime[1000];
  unsigned trigCount[N_TRIG];
  
  
  for(int irfile=0; irfile<N_FILE; irfile++)
    {
      cout << "open: " << rfile_name[irfile] << endl;
      
      TFile* rfile = new TFile(rfile_name[irfile]);
      TTree* rtree = (TTree*)rfile->Get("save");
  
      rtree->SetBranchAddress("eventID", &eventID);
      rtree->SetBranchAddress("eventType", &eventType);
      rtree->SetBranchAddress("triggerType", &triggerType);
      rtree->SetBranchAddress("nHits", &nHits);
      rtree->SetBranchAddress("boardID", boardID);
      rtree->SetBranchAddress("channelID", channelID);
      rtree->SetBranchAddress("tdcTime", tdcTime);
      rtree->SetBranchAddress("triggerTime", triggerTime);

      
      int nevt = rtree->GetEntries();
      //int nevt = 280000;
      
      cout << "total " << nevt << " events" << endl;

      for(int itrig=0; itrig<N_TRIG; itrig++)
        {
          trigCount[itrig]=0;
        }      
      
      for(int ievt=0; ievt<nevt; ievt++)
        {

          rtree->GetEntry(ievt);

          if( eventType!= 14)
            continue;
          

                  
          // ***  Initialize hit vector
          vector<short> tdc_hit[N_BOARD][N_CH];
          
          // ***  Loop over hits in an event  
          for(int ihit=0; ihit<nHits; ihit++)
            {                      
              bool valid_board = false;
              int iboard = -9999;
              
              for(iboard=0; iboard<N_BOARD; iboard++)
                {
                  if( boardID[ihit] == stoul( get<1>(info_board[iboard]), NULL, 16 ) )
                    {
                      valid_board = true;
                      break;
                    }
                }
                      
              bool valid_ch = false; 
              int ich = channelID[ihit];
                      
              //if( ich >=0 && ich <80 )
              if( ich >=0 && ich <96 )

                {
                  valid_ch = true;
                }

              // ***  Add hit timing
              if( valid_board & valid_ch )
                {
                  // tdcTime = TDC common stop timing (v1495 G1) - TDC hit timing (v1495 A-F)
                  tdc_hit[iboard][ich].push_back( tdcTime[ihit] );
                }
              else
                {
                  cout<<"Warning: Bad Data for event "<<ievt
                      <<" invalid board or channel ID: "
                      <<boardID[ihit]<<" "
                      <<channelID[ihit]<<" "
                      <<"at hit "<<ihit<<endl;
                }
            }// ***  End nhit loop

          // ***  multiplicity per event (ignore this part if you are not interested in)

          int hitMultiBoard  [N_BOARD];  // hit per board
          int hitMultiBoardStation[N_BOARD][N_STATION]; // hit per board per station

          int hitMultiAllHodo=0; // all hodoscope hit sum
          int hitMultiStationHodo[N_STATION]; // all hodoscope hit sum per station

           for(int iboard=0; iboard<N_BOARD; iboard++)
             for(int istation=0; istation<N_STATION; istation++)
               {
                 hitMultiStationHodo[istation]=0;
                 hitMultiBoard[iboard] = 0;
                 hitMultiBoardStation[iboard][istation]=0;
               }
           
          for(int iboard=0; iboard<N_BOARD; iboard++)
            {                    
              for(int ich=0; ich<N_CH; ich++)
                {
                  if( ich==64 || ich==65 )
                    continue;

                  hitMultiBoard[iboard] += tdc_hit[iboard][ich].size();

                  if( iboard<2 )
                    {
                      hitMultiAllHodo  += tdc_hit[iboard][ich].size();
                    }
                  
                  int istation = stoi(tdc[iboard][ich].station)-1;
                  
                  if( istation<N_STATION )
                    {
                      hitMultiBoardStation[iboard][istation] += tdc_hit[iboard][ich].size();

                      if( iboard<2 )
                        hitMultiStationHodo[istation] += tdc_hit[iboard][ich].size();
                    }
                }
            }
          
          // ***  Fill histo at the end of each event
          for(int itrig=0; itrig<N_TRIG; itrig++)
            {              
              if( ((triggerType >> itrig) & 1) && triggerType>0 )
                {                  
                  trigCount[itrig]++;

                  //if( itrig==5 )
                  //cout<<"event ID = "<<eventID<< "\t has NIM1 trigger"<< endl;
                  
                  
                  h_multi_all [itrig] -> Fill( hitMultiAllHodo );


                  for( int istation=0; istation<N_STATION; istation++)
                    {
                      h_multi_hst[itrig][istation] -> Fill( hitMultiStationHodo[istation] );
                    }
                  
                  for(int iboard=0; iboard<N_BOARD; iboard++)
                    {
                      h_multi [itrig][iboard] -> Fill( hitMultiBoard[iboard] );
                      
                      
                      for( int istation=0; istation<N_STATION; istation++)
                        {
                          h_multi_st  [itrig][iboard][istation] -> Fill( hitMultiBoardStation[iboard][istation] );
                        }
                      
                      for(int ich=0; ich<N_CH; ich++)
                        {
                          if( ich==64 || ich==65 )
                            continue;
                          
                          int nMulti = tdc_hit[iboard][ich].size();
                          
                          hh_multi [itrig][iboard]     -> Fill( ich, nMulti );
                          
                          if( nMulti>0 )
                            {
                              h_hit   [itrig][iboard] -> Fill( ich, nMulti);
                            }
                          
                          for(int imulti=0; imulti<nMulti; imulti++)
                            {                            
                              h_tdc    [itrig][iboard][ich] ->Fill( tdc_hit[iboard][ich][imulti]* TDC_NS_CONV );
                              hh_tdc   [itrig][iboard]      ->Fill( ich, tdc_hit[iboard][ich][imulti]* TDC_NS_CONV );
                              hh_t_m   [itrig][iboard][ich] ->Fill( nMulti, tdc_hit[iboard][ich][imulti]* TDC_NS_CONV );
                            }
                        }
                    }
                }
            }
        }// ***  End event loop
      
      
      for(int itrig=0; itrig<N_TRIG; itrig++)
        {
          cout << "Total \t" << trigCount[itrig] << "\t" << get<1>(info_trig[itrig]) << " trigger in this run " << endl;
        }      
      
      rfile->Close();
      rfile->Delete();
    }



  // ====================================================
  //
  // Get peak timing
  //
  // ====================================================
  
  
  float y_time  [N_TRIG][N_BOARD][N_CH];
  float x_ch    [N_TRIG][N_BOARD][N_CH];
  float ey_time [N_TRIG][N_BOARD][N_CH];
  float ex_ch   [N_TRIG][N_BOARD][N_CH];
  
  
  for(int itrig=0; itrig<N_TRIG; itrig++)
    for(int iboard=0; iboard<N_BOARD; iboard++)
      for(int ich=0; ich<N_CH; ich++)
        {
          if( stoi(tdc[iboard][ich].disc)<5 && iboard<4 && (itrig==5 || itrig==6 || itrig==8) )
            {
              // get_peak() is defined somewhere in the head of this code
              // get_peak() gives you the peak TDC bin for now, for cosmic
              // get_peak() should be modified to gaus fit for beam commissioning in the future
              y_time [itrig][iboard][ich]= get_peak( h_tdc[itrig][iboard][ich] );          
              ey_time[itrig][iboard][ich]= 15*TDC_NS_CONV;
              x_ch   [itrig][iboard][ich]= ich;
              ex_ch  [itrig][iboard][ich]= 1;
              
              cout << get<1>(info_trig[itrig])  <<"\t"
                   << "0x"  << get<1>(info_board[iboard]) <<"\t"
                   << "ch "              << ich <<"\t"
                   << "peak timing = "   << y_time [itrig][iboard][ich] <<"\t"
                   << "histo entries = " << h_tdc[itrig][iboard][ich]->GetEntries() <<"\t"
                   << "histo mean = "    << h_tdc[itrig][iboard][ich]->GetMean() << endl;
            }
          else
            {
              y_time [itrig][iboard][ich]= - 9999; 
              ey_time[itrig][iboard][ich]= - 9999; 
              x_ch   [itrig][iboard][ich]= - 9999; 
              ex_ch  [itrig][iboard][ich]= - 9999; 
     
            }          
        }


  // ====================================================
  //
  // Set delay timing
  //
  // ====================================================
  
  
  float new_delay [N_BOARD][N_CH];
  
  for(int iboard=0; iboard<N_BOARD; iboard++)
    {
      ifstream ifs(get<4>(info_board[iboard]));            

      for(int ich=0; ich<N_CH; ich++)
        {
          char line_delay[255];
          ifs.getline(line_delay,255);
          
          int read_ch;
          int read_delay;

          sscanf(line_delay,"%d:f%x",&read_ch,&read_delay);

          //cout << line_delay << " = " << read_ch << " :f " << read_delay << endl;
          
          new_delay[iboard][ich] = update_delay(read_delay,  y_time[5][iboard][ich]);

          // cout << "new delay value = " << read_ch << " :f " << new_delay[iboard][ich] << endl << endl;

        }
    }

  
  // ====================================================
  //
  // Timing graph
  //
  // ====================================================


  // ***  peak timing vs channel
  
  TGraphErrors* g_tdc[N_TRIG][N_BOARD];

  for(int itrig=0; itrig<N_TRIG; itrig++)
    for(int iboard=0; iboard<N_BOARD; iboard++)
      {
        string obj_id;
        obj_id = "g_trig"+to_string(itrig)+"_board"+to_string(iboard);

        string label;
        label = run_name +" " + get<1>(info_trig[itrig]) + " " + get<1>(info_board[iboard]) + "; Time (ns) ; TDC Channel";

        g_tdc[itrig][iboard] = new TGraphErrors
          (N_CH, x_ch[itrig][iboard], y_time[itrig][iboard], ex_ch[itrig][iboard], ey_time[itrig][iboard] );

        g_tdc[itrig][iboard]->Draw("ap");
        g_tdc[itrig][iboard]->SetName(obj_id.c_str());
        g_tdc[itrig][iboard]->SetTitle(label.c_str());
      }
  
    
  // ====================================================
  //
  // Drawing
  //
  // ====================================================
  

  TString pdf_tag(argv[1]);
  TString pdf_name = "pdf_files/c_tdc_"+pdf_tag+".pdf";
    
  TCanvas* c_tdc = new TCanvas("c_tdc","TDC",2400,1200);
  
  
  TString pdf_name_temp = pdf_name+"[";

  cout<<"generating pdf file: "<<pdf_name<<endl;
  
        
  c_tdc -> Print(pdf_name_temp);
  //c_tdc -> Print("c_tdc.pdf[");


  // *** #hit vs chID
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      if( itrig!=5 && itrig!= 6 && itrig!=8 )
        continue;
      
      c_tdc->Divide(2,2);

      for(int iboard=0; iboard<N_BOARD-1; iboard++)
        {        
          c_tdc->cd(iboard+1);
          h_hit[itrig][iboard] -> Sumw2(0);
          h_hit[itrig][iboard] -> Draw();
          gPad->SetLogy();
          
        }
      
      c_tdc->cd(0);
      //c_tdc -> Print("c_tdc.pdf");
      c_tdc -> Print(pdf_name);
      
      c_tdc -> Clear();
    }

  // *** #time vs chID
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      if( itrig!=5 && itrig!= 6 && itrig!=8 )
        continue;
      
      c_tdc->Divide(2,2);

      for(int iboard=0; iboard<N_BOARD-1; iboard++)
        {        
          c_tdc->cd(iboard+1);
          hh_tdc[itrig][iboard] -> Draw("colz");
          gPad->SetLogz();
        }
      
      c_tdc->cd(0);
      //c_tdc -> Print("c_tdc.pdf");
      c_tdc -> Print(pdf_name);
      c_tdc -> Clear();
    }

  /*** this part is commented out, but works ok
  
  // *** #event vs #multiplicity per board per station
  
  c_tdc->Divide(2,2);
  
  c_tdc->cd(1);
  h_multi_all[5] -> Draw();
  gPad->SetLogx();
  gPad->SetLogy();
  
  c_tdc->cd(2);
  h_multi_all[6] -> Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  c_tdc->cd(3);
  h_multi_all[7] -> Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  c_tdc->cd(4);
  h_multi_all[8] -> Draw();
  gPad->SetLogx();
  gPad->SetLogy();
                 
  c_tdc->cd(0);
  c_tdc -> Print(pdf_name);
  c_tdc -> Clear();

  // *** #event vs #multiplicity per station
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      if( itrig!=5 && itrig!= 6 && itrig!=8 )
        continue;
      
      c_tdc->Divide(2,2);
      
      for(int istation=0; istation<N_STATION; istation++)
        {
          c_tdc->cd(istation+1);
          h_multi_hst[itrig][istation] -> Draw();
           gPad->SetLogx();
           gPad->SetLogy();
        }
      
      c_tdc->cd(0);
      c_tdc -> Print(pdf_name);
      c_tdc -> Clear();
      
    }  

  // *** #event vs #multiplicity per board
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      if( itrig!=5 && itrig!= 6 && itrig!=8 )
        continue;
      
      c_tdc->Divide(2,2);

      for(int iboard=0; iboard<N_BOARD-1; iboard++)
        {        
          c_tdc->cd(iboard+1);
          h_multi[itrig][iboard] -> Draw();
          gPad->SetLogx();
          gPad->SetLogy();
        }
      
      c_tdc->cd(0);
      c_tdc -> Print(pdf_name);
      c_tdc -> Clear();
    }

  // *** #event vs #multiplicity per board per station
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      if( itrig!=5 && itrig!= 6 && itrig!=8 )
        continue;
      
      for(int iboard=0; iboard<N_BOARD-1; iboard++)
        {
          c_tdc->Divide(2,2);

          for(int istation=0; istation<N_STATION; istation++)
            {
              c_tdc->cd(istation+1);
              h_multi_st[itrig][iboard][istation] -> Draw();
              gPad->SetLogx();
              gPad->SetLogy();
            }
      
          c_tdc->cd(0);
          c_tdc -> Print(pdf_name);
          c_tdc -> Clear();
        }
    }
  
  // *** #multiplicity vs chID
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      if( itrig!=5 && itrig!= 6 && itrig!=8 )
        continue;
      
      c_tdc->Divide(2,2);

      for(int iboard=0; iboard<N_BOARD-1; iboard++)
        {        
          c_tdc->cd(iboard+1);
          hh_multi[itrig][iboard] -> Draw("colz");
          gPad->SetLogz();
        }
      
      c_tdc->cd(0);
      c_tdc -> Print(pdf_name);
      c_tdc -> Clear();
    }

  // *** peak timing vs chID
  for(int itrig=0; itrig<N_TRIG; itrig++)
    {
      if( itrig!=5 && itrig!= 6 && itrig!=8 )
        continue;
      
      c_tdc->Divide(2,2);

      for(int iboard=0; iboard<N_BOARD-1; iboard++)
        {        
          c_tdc->cd(iboard+1);
          g_tdc[itrig][iboard  ]->Draw("ap");
          g_tdc[itrig][iboard  ]->GetYaxis()->SetRangeUser(TIME_MIN,TIME_MAX);
          g_tdc[itrig][iboard  ]->GetXaxis()->SetRangeUser(0,N_CH);
          g_tdc[itrig][iboard  ]->Draw("ap");
        }
      
      c_tdc->cd(0);
      c_tdc -> Print(pdf_name);
      c_tdc -> Clear();
    }
  **/

  pdf_name_temp = pdf_name+"]";
        
  c_tdc -> Print(pdf_name_temp);
  //c_tdc -> Print("c_tdc.pdf]");
  
  cout << "done" << endl;

    
}
