/*
================================================================================
Function: makeClass()
    Reads the ROOT files exist in the playlist
    Outputs Truth_Branch.C and Truth_Branch.h files
    Truth_Branch.h contains member variables for the ROOT File
    Use Truth.C to analyze ALL Truth events

    Usage:
        > .L makeClass.cpp
        > makeClass()
    
================================================================================
*/

void makeClass(){

    string playlist = "../../Input/Playlists/pl_MC_Merged.dat";
    TChain fChain("Truth") ;

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

    fChain.MakeClass("Truth_Branch");

}
