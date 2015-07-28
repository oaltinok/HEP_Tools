/*
================================================================================
Function: makeClass()
    Reads the ROOT files exist in the playlist
    Outputs MC_Sample.C and MC_Sample.h files
    MC_Sample.h contains member variables for the ROOT File
    Copy these member variables and branch addresses to Corresponding Class
    
    Usage:
        > .L makeClass.cpp
        > makeClass()
================================================================================
*/

void makeClass(){

    string playlist = "../../Input/Playlists/pl_MC_Merged.dat";
    TChain fChain("CCProtonPi0") ;

    ifstream input_pl( playlist.c_str() );
    string filename;

    if( !input_pl.is_open() ){
        cerr<<"Cannot open Playlist File!"<<endl;
        exit(1);
    }else{
        cout<<"Playlist: "<<playlist.c_str()<<endl;
    }

   while (input_pl) {
     input_pl>>filename;
     
     if (!input_pl) break;
    
     if (filename[0] != '/') break;
    
     fChain.Add( filename.c_str() );
     cout<<" Added "<<filename.c_str()<<endl;
   }

    fChain.MakeClass("MC_Sample");

}
