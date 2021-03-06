/*
================================================================================
Function: makeClass()
    Creates a Class for Meta Tree
================================================================================
*/

void makeClass(){

    string playlist = "../../Input/Playlists/pl_MC_Merged.dat";
    TChain fChain("Meta") ;

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

    fChain.MakeClass("Meta_Class");

}
