/* Build charge, pedestal, and average waveform histograms
   from oscilloscope data in .hdf5 format
   * * *
   * input params:
   * argv[1]: text file listing .hdf5 for analysis
   * argv[2]: name of output root file to write histograms to
   * argv[3]: length of the pedestal window
   * argv[4]: length of the integration window
   * argv[5]: scope channel for analysis
   * * *
*/

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>

// Root Libraries
#include <TROOT.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TPaveStats.h>
#include <TLatex.h>

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

// Initialize Datacluster
DataCluster * Init_Data(DataSet *dataset);

// For each trace, read in the dataset
int Read_Trace(DataCluster *datacluster, unsigned long trace_index);

// Return the voltage of the time-bin
double getVoltage(DataCluster *datacluster, double pedestal, double dy, int i);

// Return the integrated charge over the signal window
double getCharge(float window_width, float signal_window, DataCluster *datacluster,
                 double pedestal, double dy, double dx);

// Return the integral of the charge histogram above 5pc
void chargeIntegral(TH1F* charges_signal);

// Beautify charge plot
void prettyPlot(TH1F* h);

const int termination_ohms = 50;

int main (int argc, char* argv[])
{
    /* User inputs */
    float window_width = atoi(argv[3]); // pedestal window
    float signal_window = atoi(argv[4]); // signal window
    string channel = argv[5]; // analysis channel

    std::vector<float> waveform_voltage;

    try{

        // Attribute Variables
        Attribute horiz_interval;
        Attribute vertical_gain;

        // Variables for time binning and veritcal resolution
        double dx,dy;
        float trace_count = 0.0;

        // ROOT histograms
        TH1F *pedestals = new TH1F("Pedestal","",5000,0.0,2);
        TH1F *charges_signal = new TH1F("Charge","",1250,-5.0,600);

        string filename;
        H5File file;
        DataSet dataset;

        // Load in the datafiles from a txt file
        ifstream ifs (argv[1] , ifstream::in);
        ifs >> filename;

        while (ifs.good()){

            // Open HDF5 File, then HDF5 DataSet and Read in Attributes
            file.openFile(filename, H5F_ACC_RDONLY);
            ifs >> filename;
            dataset = file.openDataSet(channel);

            horiz_interval = dataset.openAttribute("horiz_interval");
            vertical_gain = dataset.openAttribute("vertical_gain");

            horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
            vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);

            DataCluster * datacluster = Init_Data(&dataset);
            unsigned long window_length = datacluster->trace_length;

            cout << "Analyzing file:           "  << filename << endl;

            // Print the waveform parameters to screen
            cout << "Waveform parameters for   " << channel << endl;
            cout << " * * * * * * * * * * * * * * * * * * * " << endl;
            cout << "Horizontal resolution:    " << dx << " ns" << endl;
            cout << "Verticle resolution:      " << dy << " V" << endl;
            cout << "Trace length:             " << window_length*dx*1e9 << " ns" << endl;
            cout << "Pedestal window:          " << "0 -" << window_width*dx*1e9 << " ns"<< endl;
            cout << "Integration window:       " << window_width*dx*1e9 << " - "
                 << (window_width + signal_window)*dx*1e9 << " ns" << endl;

            // Loop over every trace
            for (unsigned int j = 0; j < datacluster->n_traces; j++){
                cout << "Analyzing trace " << j+1 << " out of "
                     <<datacluster->n_traces<< " total traces " <<'\r';
                cout.flush();

                // Read in the data
                Read_Trace(datacluster,j);

                double pedestal = TMath::Mean (window_width, datacluster->data_out)*dy; // baseline
                unsigned int window_length = datacluster->trace_length; // trace length

                if(j == 0){ // resize once
                    waveform_voltage.resize(window_length);
                }

                double charge = getCharge(window_width, signal_window, datacluster, pedestal, dy, dx);

                for(unsigned int i = 0; i < window_length; i++){
                    waveform_voltage[i] += getVoltage(datacluster, pedestal, dy, i);
                }

                charges_signal->Fill(charge);
                pedestals->Fill(pedestal);
                trace_count++;
            }

            file.close();
        }

        ifs.close(); // close file stream

        // Build the average waveform
        TH1F* avg_wfm = new TH1F("avg_wfvm","",waveform_voltage.size(),0,waveform_voltage.size()*dx*1e9);
        for(unsigned int i = 0; i < waveform_voltage.size(); i++){
            avg_wfm->SetBinContent(i, waveform_voltage[i]/trace_count);
        }

        chargeIntegral(charges_signal);
        prettyPlot(charges_signal);

        // Output Histograms to File
        if (argc > 2){
            TFile f(argv[2],"new");
            pedestals->Write();
            charges_signal->Write();
            avg_wfm->Write();
        }

        // clean-up
        delete pedestals;
        delete charges_signal;
        delete avg_wfm;

    } // end of try block

    // catch failure caused by the H5File operations
    catch(FileIException error){
        error.printError();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error){
        error.printError();
        return -1;
    }

    // catch failure caused by the DataSpace operations
    catch(DataSpaceIException error){
        error.printError();
        return -1;
    }

    // catch failure caused by the DataSpace operations
    catch(DataTypeIException error){
        error.printError();
        return -1;
    }

    return 0; // successfully terminated
}

/* Datacluster method */
DataCluster * Init_Data(DataSet *dataset){

    DataCluster * datacluster = new DataCluster[1];

    /* Get dataspace of the dataset. */
    datacluster->dataset = dataset;
    datacluster->dataspace = datacluster->dataset->getSpace();

    /* Get the dimension size of each dimension
       in the dataspace and display them. */
    hsize_t dims_out[2];
    datacluster->dataspace.getSimpleExtentDims(dims_out, NULL);
    datacluster->trace_length = (unsigned long)(dims_out[1]);
    datacluster->n_traces = (unsigned long)(dims_out[0]);

    // Data Buffer
    datacluster->data_out = new char[datacluster->trace_length]; // Data is size char
    for (unsigned long i = 0; i < datacluster->trace_length; i++) datacluster->data_out[i]= 0;

    /* Define hyperslab in the dataset. */
    datacluster->offset[0] = 0;
    datacluster->offset[1] = 0;
    datacluster->count[0] = 1;
    datacluster->count[1] = datacluster->trace_length;
    datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );

    /* Define the memory dataspace. */
    hsize_t dimsm[2]; /* memory space dimensions */
    dimsm[0] = dims_out[0];
    dimsm[1] = dims_out[1];
    datacluster->memspace = DataSpace( RANK_OUT, dimsm );

    /* Define memory hyperslab. */
    datacluster->offset_out[0] = 0;
    datacluster->offset_out[1] = 0;
    datacluster->count_out[0] = 1;
    datacluster->count_out[1] = datacluster->trace_length;
    datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );

    return datacluster;
}

/* Method: Read_Trace(DataCluster *datacluster, unsigned long trace_index)
 * Updates a DataCluster datacluster so that its buffer contains trace number trace_index */
int Read_Trace(DataCluster *datacluster, unsigned long trace_index){
    datacluster->offset[0]= (hsize_t)trace_index;
    datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );
    datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );
    datacluster->dataset->read( datacluster->data_out, PredType::NATIVE_CHAR, datacluster->memspace, datacluster->dataspace );
    return 0;
}

// Returns voltage of the time-bin
double getVoltage(DataCluster *datacluster, double pedestal, double dy, int i){
    float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
    return voltage;
}

// Returns charge (pC) integrated over the signal window
double getCharge(float window_width, float signal_window, DataCluster *datacluster,
                 double pedestal, double dy, double dx){
    double charge = 0.0;
    for(int i = window_width; i < window_width + signal_window; i++){
        float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
        charge+=(voltage*((-1000.0*dx*1e9)/termination_ohms)); // in pC
    }

    return charge;
}

// Charge-weighted integral of the charge histogram
void chargeIntegral(TH1F* charges_signal){

    TAxis *axis = charges_signal->GetXaxis();

    // Integration range
    int bmin = axis->FindBin(5.0);
    int bmax = axis->FindBin(600.0);

    double x,y,integral = 0.0;

    // Charge-weighted integral of charge histogram
    for(int i = bmin; i < bmax; i++){
        x = axis->GetBinCenter(i);
        y = charges_signal->GetBinContent(i);
        integral += x*y;
    }

    // Print the integral of the charge distribtuion above some noise threshold
    cout << "Weighted integral above " << 5.0 << "pC is " << integral << endl;

    TCanvas *c = new TCanvas;
    charges_signal->Draw();
    c->Update();

}

// Beautify charge plot
void prettyPlot(TH1F* h){

    TAxis *xaxis = h->GetXaxis();
    TAxis *yaxis = h->GetYaxis();
    xaxis->SetTitleFont(132);
    xaxis->SetLabelFont(132);
    xaxis->SetTitle("Charge (pC)");
    yaxis->SetTitleFont(132);
    yaxis->SetLabelFont(132);
    yaxis->SetTitle("Traces");
    h->SetLineColor(kBlack);

}

