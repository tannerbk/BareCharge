#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

// Root Libraries
#include <TFile.h>
#include <TMath.h>
#include <TRandom2.h>
#include <TMathBase.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

// HDF5 Library
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

const int RANK_OUT = 2; // Data Rank

struct DataCluster{
  DataSet *dataset; // Dataset pointer
  DataSpace dataspace; // DataSet's DataSpace
  DataSpace memspace; // MemSpace Object for Data Extraction
  hsize_t offset[RANK_OUT]; // Data Extraction Parameters...
  hsize_t count[RANK_OUT];
  hsize_t offset_out[RANK_OUT];
  hsize_t count_out[RANK_OUT];
  unsigned long trace_length; // Length of a Scope Trace 
  unsigned long n_traces; // Number of traces in DataSet
  char * data_out; // Pointer to Data Buffer
};

typedef struct DataCluster DataCluster;

/* DataCluster Methods */

DataCluster * Init_Data(DataSet *dataset);
int Read_Trace(DataCluster *datacluster, unsigned long trace_index);

int main (int argc, char* argv[])
{
  unsigned long window_width = atoi(argv[3]);
  const int termination_ohms = 50; 

  try
    {

      // Attribute Variables
      Attribute horiz_interval;
      Attribute vertical_gain;
      Attribute horiz_offset;
      Attribute vertical_offset;
      Attribute max_value;
      Attribute min_value;
      
      // Attribute Value Variables
      double dx,dy,xoffset,yoffset;
      double Vmin, Vmax;

      //Define ROOT histograms here, the variable names are confusing right now, look at the descriptions 
      TH1F *pedestals = new TH1F("Pedestal","",5000,0.5,2); 
      TH1F *variances = new TH1F("Variance over Pedestal Window","",5000,-0.1,1.0);
      TH1F *charges_signal = new TH1F("Charge PMT","",100000,-150,300);   
      TH1F *average_waveform = new TH1F("Average Waveform","", 5002, 0, 500); 
	          
      unsigned long i; 
      unsigned long j; 
      unsigned int window_count = 0.0; // counts number of traces 
      string filename; 
      H5File file; 
      DataSet dataset;

      const unsigned long length_trace = 10002; // length of trace  
      float waveform_voltage[length_trace] = {0};

      ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
      ifs >> filename;

      while (ifs.good()){ // While Stream Is Open, Analyze New Files.

	cout<<filename<<endl;
	file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File 
	ifs >> filename;
	dataset = file.openDataSet("channel2"); // Open HDF5 DataSet, Channel 2

	// Reading in Attributes
	horiz_interval = dataset.openAttribute("horiz_interval");
	vertical_gain = dataset.openAttribute("vertical_gain");
	horiz_offset = dataset.openAttribute("horiz_offset");
	vertical_offset = dataset.openAttribute("vertical_offset");
	max_value = dataset.openAttribute("max_value");
	min_value = dataset.openAttribute("min_value");

	horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
	vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);
	horiz_offset.read(PredType::NATIVE_DOUBLE, &xoffset);
	vertical_offset.read(PredType::NATIVE_DOUBLE, &yoffset);
	max_value.read(PredType::NATIVE_DOUBLE, &Vmax);
	min_value.read(PredType::NATIVE_DOUBLE, &Vmin);
	  
	DataCluster * datacluster = Init_Data(&dataset); //initialize datacluser
	
	unsigned long window_length = datacluster->trace_length; 
	cout << " The time bin width is " << dx  << ", the trace length is " << window_length << ", the verticle resolution is " << dy << "." << endl;
		 
	double ncharge; 
	float voltage = 0.0; 
	float signal_voltage = 0.0; 
	float window_voltage = 0.0; 
	float variance; 
	float pedestal = 0.0; 

	for (j = 0; j < datacluster->n_traces; j++){ // loop over all the traces 

	  cout << "Analyzing trace " << j+1 << " out of " 
	       <<datacluster->n_traces<< " total traces " <<'\r';
	  cout.flush(); 

	  Read_Trace(datacluster,j); 	    	    
	  pedestal = TMath::Mean (window_width, &datacluster->data_out[window_width])*dy; // calculates pedestal 
	  
	  ncharge = 0.0;  
	  variance = 0.0; 	   
	    
	  for(i=window_width; i < window_length; i++){ // calculates variance in second half of the window 
	    voltage = ((float)datacluster->data_out[i]*dy-pedestal);
	    variance = variance + voltage * voltage; 
	  }
	    
	  // no variance cut, ie not cutting dark pulses from pedestal window. variance cut could be put in here if there is a concern about dark pulses.  

	  for(i = 900; i < window_width + 1200 ; i++){ // signal window, hardcoded
	    signal_voltage = ((float)datacluster->data_out[i]*dy-pedestal);
	    ncharge = ncharge+(signal_voltage*((-1000.0*dx*1e9)/termination_ohms)); // charge 
	  }  

	  for(i=0; i < window_length; i++){ // entire window, build array for average waveform 
	    window_voltage = ((float)datacluster->data_out[i]*dy-pedestal);  
	    waveform_voltage[i] = waveform_voltage[i] + window_voltage; 
	  }	   

	  window_count++; // number of traces 
	  
	  // fill histograms 
	  charges_signal->Fill(ncharge); 
	  variances->Fill(variance); 
	  pedestals->Fill(pedestal);
	}
	
	file.close(); // close file 
      }
      ifs.close(); // close file stream 
    
      // finds averge waveform 
      Float_t avg_waveform_voltage[length_trace];  
      for(i=0; i<length_trace; i++){ 
	avg_waveform_voltage[i] = waveform_voltage[i]/window_count;
	average_waveform->SetBinContent(i+1, avg_waveform_voltage[i]);  
      }

      // Output Histograms to File
      if (argc > 2){
	TFile f(argv[2],"new");
	pedestals->Write();
	variances->Write();
	charges_signal->Write();
	average_waveform->Write();

      }

      // clean-up 
      delete pedestals; 
      delete variances; 
      delete charges_signal;   
      delete average_waveform; 

    } // end of try block
  
  // catch failure caused by the H5File operations
  catch( FileIException error )
    {
      error.printError();
      return -1;
    }

  // catch failure caused by the DataSet operations
  catch( DataSetIException error )
    {
      error.printError();
      return -1;
    }

  // catch failure caused by the DataSpace operations
  catch( DataSpaceIException error )
    {
      error.printError();
      return -1;
    }

  // catch failure caused by the DataSpace operations
  catch( DataTypeIException error )
    {
      error.printError();
      return -1;
    }

  return 0; // successfully terminated
}

DataCluster * Init_Data(DataSet *dataset){

  DataCluster * datacluster = new DataCluster[1];
  
  /*
   * Get dataspace of the dataset.
   */
  datacluster->dataset = dataset;
  datacluster->dataspace = datacluster->dataset->getSpace();
  
  /*
   * Get the dimension size of each dimension in the dataspace and
   * display them.
   */
  
  hsize_t dims_out[2];
  datacluster->dataspace.getSimpleExtentDims( dims_out, NULL);
  datacluster->trace_length = (unsigned long)(dims_out[1]);
  datacluster->n_traces = (unsigned long)(dims_out[0]);
      
  // cout << "Reading " << datacluster->n_traces << " traces of length " << datacluster->trace_length << "..." << endl;
  
  // Data Buffer
  datacluster->data_out = new char[datacluster->trace_length]; // Scope data is size char
  for (unsigned long i = 0; i < datacluster->trace_length; i++) datacluster->data_out[i]= 0;

  /*
   * Define hyperslab in the dataset.
   */

  datacluster->offset[0] = 0;
  datacluster->offset[1] = 0;
  datacluster->count[0] = 1;
  datacluster->count[1] = datacluster->trace_length;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );
  
  /*
   * Define the memory dataspace.
   */

  hsize_t dimsm[2]; /* memory space dimensions */
  dimsm[0] = dims_out[0];
  dimsm[1] = dims_out[1];
  datacluster->memspace = DataSpace( RANK_OUT, dimsm );

  /*
   * Define memory hyperslab.
   */

  datacluster->offset_out[0] = 0;
  datacluster->offset_out[1] = 0;
  datacluster->count_out[0] = 1;
  datacluster->count_out[1] = datacluster->trace_length;
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );
  
  return datacluster;
}

/* Method: Read_Trace(DataCluster *datacluster, unsigned long trace_index)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Updates a DataCluster datacluster so that its buffer contains trace number trace_index
 */

int Read_Trace(DataCluster *datacluster, unsigned long trace_index){
  datacluster->offset[0]= (hsize_t)trace_index;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );
  datacluster->dataset->read( datacluster->data_out, PredType::NATIVE_CHAR, datacluster->memspace, datacluster->dataspace );
  return 0; // No protection...
}
