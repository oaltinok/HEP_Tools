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
    
    Last Revision: 2014_01_06
================================================================================
*/

void makeClass(){

    string playlist = "Playlists/pl_DST_01.dat";
    TChain fChain("minerva") ;

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

    fChain.MakeClass("PC_DST_Sample");

}
