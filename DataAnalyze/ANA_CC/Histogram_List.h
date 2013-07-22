    //------------------------------------------------------------------------
    //  Truth vs Reco Comparisons
    //------------------------------------------------------------------------

    // Incoming Neutrino Energy

    TH2F* Ev_reco_Ev_true_t0 = new TH2F("Ev_reco_Ev_true_t0","E_{#nu} True vs E_{#nu} Reco for 0 Pions",
        NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_Ev, MIN_Ev, MAX_Ev );
    Ev_reco_Ev_true_t0->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_reco_Ev_true_t0->GetYaxis()->SetTitle("True Neutrino Energy [GeV]");

    TH2F* Ev_reco_Ev_true_t1 = new TH2F("Ev_reco_Ev_true_t1","E_{#nu} True vs E_{#nu} Reco for 1 Pions",
        NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_Ev, MIN_Ev, MAX_Ev );
    Ev_reco_Ev_true_t1->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_reco_Ev_true_t1->GetYaxis()->SetTitle("True Neutrino Energy [GeV]");

    TH2F* Ev_reco_Ev_true_t2 = new TH2F("Ev_reco_Ev_true_t2","E_{#nu} True vs E_{#nu} Reco for 2 Pions",
        NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_Ev, MIN_Ev, MAX_Ev );
    Ev_reco_Ev_true_t2->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_reco_Ev_true_t2->GetYaxis()->SetTitle("True Neutrino Energy [GeV]");

    TH2F* Ev_reco_Ev_true_t3 = new TH2F("Ev_reco_Ev_true_t3","E_{#nu} True vs E_{#nu} Reco for 3 Pions",
        NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_Ev, MIN_Ev, MAX_Ev );
    Ev_reco_Ev_true_t3->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_reco_Ev_true_t3->GetYaxis()->SetTitle("True Neutrino Energy [GeV]");

    TH2F* Ev_reco_Ev_true_t4 = new TH2F("Ev_reco_Ev_true_t4","E_{#nu} True vs E_{#nu} Reco for 4+ Pions",
        NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_Ev, MIN_Ev, MAX_Ev );
    Ev_reco_Ev_true_t4->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_reco_Ev_true_t4->GetYaxis()->SetTitle("True Neutrino Energy [GeV]");


    // Q-Square

    TH2F* q_reco_q_true_t0 = new TH2F("q_reco_q_true_t0","Q^{2} True vs Q^{2} Reco for 0 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_reco_q_true_t0->GetXaxis()->SetTitle("Q^{2} Reco [GeV]");
    q_reco_q_true_t0->GetYaxis()->SetTitle("Q^{2} True [GeV]");

    TH2F* q_reco_q_true_t1 = new TH2F("q_reco_q_true_t1","Q^{2} True vs Q^{2} Reco for 1 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_reco_q_true_t1->GetXaxis()->SetTitle("Q^{2} Reco [GeV]");
    q_reco_q_true_t1->GetYaxis()->SetTitle("Q^{2} True [GeV]");

    TH2F* q_reco_q_true_t2 = new TH2F("q_reco_q_true_t2","Q^{2} True vs Q^{2} Reco for 2 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_reco_q_true_t2->GetXaxis()->SetTitle("Q^{2} Reco [GeV]");
    q_reco_q_true_t2->GetYaxis()->SetTitle("Q^{2} True [GeV]");

    TH2F* q_reco_q_true_t3 = new TH2F("q_reco_q_true_t3","Q^{2} True vs Q^{2} Reco for 3 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_reco_q_true_t3->GetXaxis()->SetTitle("Q^{2} Reco [GeV]");
    q_reco_q_true_t3->GetYaxis()->SetTitle("Q^{2} True [GeV]");

    TH2F* q_reco_q_true_t4 = new TH2F("q_reco_q_true_t4","Q^{2} True vs Q^{2} Reco for 4+ Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_reco_q_true_t4->GetXaxis()->SetTitle("Q^{2} Reco [GeV]");
    q_reco_q_true_t4->GetYaxis()->SetTitle("Q^{2} True [GeV]");


    // W

    TH2F* w_reco_w_true_t0 = new TH2F("w_reco_w_true_t0","W True vs W Reco for 0 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    w_reco_w_true_t0->GetXaxis()->SetTitle("W Reco [GeV]");
    w_reco_w_true_t0->GetYaxis()->SetTitle("W True [GeV]");

    TH2F* w_reco_w_true_t1 = new TH2F("w_reco_w_true_t1","W True vs W Reco for 1 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    w_reco_w_true_t1->GetXaxis()->SetTitle("W Reco [GeV]");
    w_reco_w_true_t1->GetYaxis()->SetTitle("W True [GeV]");

    TH2F* w_reco_w_true_t2 = new TH2F("w_reco_w_true_t2","W True vs W Reco for 2 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    w_reco_w_true_t2->GetXaxis()->SetTitle("W Reco [GeV]");
    w_reco_w_true_t2->GetYaxis()->SetTitle("W True [GeV]");

    TH2F* w_reco_w_true_t3 = new TH2F("w_reco_w_true_t3","W True vs W Reco for 3 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    w_reco_w_true_t3->GetXaxis()->SetTitle("W Reco [GeV]");
    w_reco_w_true_t3->GetYaxis()->SetTitle("W True [GeV]");

    TH2F* w_reco_w_true_t4 = new TH2F("w_reco_w_true_t4","W True vs W Reco for 4+ Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    w_reco_w_true_t4->GetXaxis()->SetTitle("W Reco [GeV]");
    w_reco_w_true_t4->GetYaxis()->SetTitle("W True [GeV]");


    //------------------------------------------------------------------------
    //  1D Histograms
    //------------------------------------------------------------------------

    // Incoming Neutrino Energy Reco

    TH1F *Ev_reco_t0 = new TH1F("Ev_reco_t0","Incoming Neutrino Energy for 0 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_reco_t0->GetXaxis()->SetTitle("Reconstructed Neutrino Energy E_{#nu} [GeV]");
    Ev_reco_t0->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_reco_t1 = new TH1F("Ev_reco_t1","Incoming Neutrino Energy for 1 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_reco_t1->GetXaxis()->SetTitle("Reconstructed Neutrino Energy E_{#nu} [GeV]");
    Ev_reco_t1->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_reco_t2 = new TH1F("Ev_reco_t2","Incoming Neutrino Energy for 2 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_reco_t2->GetXaxis()->SetTitle("Reconstructed Neutrino Energy E_{#nu} [GeV]");
    Ev_reco_t2->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_reco_t3 = new TH1F("Ev_reco_t3","Incoming Neutrino Energy for 3 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_reco_t3->GetXaxis()->SetTitle("Reconstructed Neutrino Energy E_{#nu} [GeV]");
    Ev_reco_t3->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_reco_t4 = new TH1F("Ev_reco_t4","Incoming Neutrino Energy for 4+ Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_reco_t4->GetXaxis()->SetTitle("Reconstructed Neutrino Energy E_{#nu} [GeV]");
    Ev_reco_t4->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    // Incoming Neutrino Energy True

    TH1F *Ev_true_t0 = new TH1F("Ev_true_t0","Incoming Neutrino Energy for 0 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_true_t0->GetXaxis()->SetTitle("True Neutrino Energy E_{#nu} [GeV]");
    Ev_true_t0->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_true_t1 = new TH1F("Ev_true_t1","Incoming Neutrino Energy for 1 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_true_t1->GetXaxis()->SetTitle("True Neutrino Energy E_{#nu} [GeV]");
    Ev_true_t1->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_true_t2 = new TH1F("Ev_true_t2","Incoming Neutrino Energy for 2 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_true_t2->GetXaxis()->SetTitle("True Neutrino Energy E_{#nu} [GeV]");
    Ev_true_t2->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_true_t3 = new TH1F("Ev_true_t3","Incoming Neutrino Energy for 3 Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_true_t3->GetXaxis()->SetTitle("True Neutrino Energy E_{#nu} [GeV]");
    Ev_true_t3->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );

    TH1F *Ev_true_t4 = new TH1F("Ev_true_t4","Incoming Neutrino Energy for 4+ Pions",NBINS_Ev, MIN_Ev, MAX_Ev);
    Ev_true_t4->GetXaxis()->SetTitle("True Neutrino Energy E_{#nu} [GeV]");
    Ev_true_t4->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_Ev) );


    
    // Reconstructed W

    TH1F *w_reco_t0 = new TH1F("w_reco_t0","Reconstructed W for 0 Pions",NBINS_W, MIN_W, MAX_W);
    w_reco_t0->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w_reco_t0->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_reco_t1 = new TH1F("w_reco_t1","Reconstructed W for 1 Pions",NBINS_W, MIN_W, MAX_W);
    w_reco_t1->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w_reco_t1->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_reco_t2 = new TH1F("w_reco_t2","Reconstructed W for 2 Pions",NBINS_W, MIN_W, MAX_W);
    w_reco_t2->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w_reco_t2->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_reco_t3 = new TH1F("w_reco_t3","Reconstructed W for 3 Pions",NBINS_W, MIN_W, MAX_W);
    w_reco_t3->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w_reco_t3->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_reco_t4 = new TH1F("w_reco_t4","Reconstructed W for 4+ Pions",NBINS_W, MIN_W, MAX_W);
    w_reco_t4->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w_reco_t4->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );


    // True W
    
    TH1F *w_true_t0 = new TH1F("w_true_t0","True W for 0 Pions",NBINS_W, MIN_W, MAX_W);
    w_true_t0->GetXaxis()->SetTitle("True W [GeV]");
    w_true_t0->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_true_t1 = new TH1F("w_true_t1","True W for 1 Pions",NBINS_W, MIN_W, MAX_W);
    w_true_t1->GetXaxis()->SetTitle("True W [GeV]");
    w_true_t1->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_true_t2 = new TH1F("w_true_t2","True W for 2 Pions",NBINS_W, MIN_W, MAX_W);
    w_true_t2->GetXaxis()->SetTitle("True W [GeV]");
    w_true_t2->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_true_t3 = new TH1F("w_true_t3","True W for 3 Pions",NBINS_W, MIN_W, MAX_W);
    w_true_t3->GetXaxis()->SetTitle("True W [GeV]");
    w_true_t3->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *w_true_t4 = new TH1F("w_true_t4","True W for 4+ Pions",NBINS_W, MIN_W, MAX_W);
    w_true_t4->GetXaxis()->SetTitle("True W [GeV]");
    w_true_t4->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    // Reconstructed Q-Square

    TH1F *q_reco_t0 = new TH1F("q_reco_t0","Reconstructed Q^{2} for 0 Pions",NBINS_W, MIN_W, MAX_W);
    q_reco_t0->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV]");
    q_reco_t0->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_reco_t1 = new TH1F("q_reco_t1","Reconstructed Q^{2} for 1 Pions",NBINS_W, MIN_W, MAX_W);
    q_reco_t1->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV]");
    q_reco_t1->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_reco_t2 = new TH1F("q_reco_t2","Reconstructed Q^{2} for 2 Pions",NBINS_W, MIN_W, MAX_W);
    q_reco_t2->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV]");
    q_reco_t2->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_reco_t3 = new TH1F("q_reco_t3","Reconstructed Q^{2} for 3 Pions",NBINS_W, MIN_W, MAX_W);
    q_reco_t3->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV]");
    q_reco_t3->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_reco_t4 = new TH1F("q_reco_t4","Reconstructed Q^{2} for 4+ Pions",NBINS_W, MIN_W, MAX_W);
    q_reco_t4->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV]");
    q_reco_t4->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    // True Q-Square

    TH1F *q_true_t0 = new TH1F("q_true_t0","TrueQ^{2} for 0 Pions",NBINS_W, MIN_W, MAX_W);
    q_true_t0->GetXaxis()->SetTitle("True Q^{2} [GeV]");
    q_true_t0->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_true_t1 = new TH1F("q_true_t1","TrueQ^{2} for 1 Pions",NBINS_W, MIN_W, MAX_W);
    q_true_t1->GetXaxis()->SetTitle("True Q^{2} [GeV]");
    q_true_t1->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_true_t2 = new TH1F("q_true_t2","TrueQ^{2} for 2 Pions",NBINS_W, MIN_W, MAX_W);
    q_true_t2->GetXaxis()->SetTitle("True Q^{2} [GeV]");
    q_true_t2->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_true_t3 = new TH1F("q_true_t3","TrueQ^{2} for 3 Pions",NBINS_W, MIN_W, MAX_W);
    q_true_t3->GetXaxis()->SetTitle("True Q^{2} [GeV]");
    q_true_t3->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );

    TH1F *q_true_t4 = new TH1F("q_true_t4","True Q^{2} for 4+ Pions",NBINS_W, MIN_W, MAX_W);
    q_true_t4->GetXaxis()->SetTitle("True Q^{2} [GeV]");
    q_true_t4->GetYaxis()->SetTitle( Form("Candidates / %3.1f ",WIDTH_W) );


    //------------------------------------------------------------------------
    //  2D Comparison Histograms
    //------------------------------------------------------------------------


    // Neutrino Energy vs W

    TH2F* Ev_w_t0 = new TH2F("Ev_w_t0","E_{#nu} Reco vs W for 0 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_w_t0->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_w_t0->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* Ev_w_t1 = new TH2F("Ev_w_t1","E_{#nu} Reco vs W for 1 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_w_t1->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_w_t1->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* Ev_w_t2 = new TH2F("Ev_w_t2","E_{#nu} Reco vs W for 2 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_w_t2->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_w_t2->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* Ev_w_t3 = new TH2F("Ev_w_t3","E_{#nu} Reco vs W for 3 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_w_t3->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_w_t3->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* Ev_w_t4 = new TH2F("Ev_w_t4","E_{#nu} Reco vs W for 4+ Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_w_t4->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_w_t4->GetYaxis()->SetTitle("Reconstructed W [GeV]");



    // Neutrino Energy vs Q-Squre

    TH2F* Ev_q_t0 = new TH2F("Ev_q_t0","E_{#nu} Reco vs Q^{2} for 0 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_q_t0->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_q_t0->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV]");

    TH2F* Ev_q_t1 = new TH2F("Ev_q_t1","E_{#nu} Reco vs Q^{2} for 1 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_q_t1->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_q_t1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV]");

    TH2F* Ev_q_t2 = new TH2F("Ev_q_t2","E_{#nu} Reco vs Q^{2} for 2 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_q_t2->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_q_t2->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV]");

    TH2F* Ev_q_t3 = new TH2F("Ev_q_t3","E_{#nu} Reco vs Q^{2} for 3 Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_q_t3->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_q_t3->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV]");

    TH2F* Ev_q_t4 = new TH2F("Ev_q_t4","E_{#nu} Reco vs Q^{2} for 4+ Pions",NBINS_Ev, MIN_Ev, MAX_Ev, NBINS_W, MIN_W, MAX_W);
    Ev_q_t4->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    Ev_q_t4->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV]");



    // Q^2 vs W

    TH2F* q_w_t0 = new TH2F("q_w_t0","Q^2 vs W for 0 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_w_t0->GetXaxis()->SetTitle("Reconstructed Q^2 [GeV]");
    q_w_t0->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* q_w_t1 = new TH2F("q_w_t1","Q^2 vs W for 1 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_w_t1->GetXaxis()->SetTitle("Reconstructed Q^2 [GeV]");
    q_w_t1->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* q_w_t2 = new TH2F("q_w_t2","Q^2 vs W for 2 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_w_t2->GetXaxis()->SetTitle("Reconstructed Q^2 [GeV]");
    q_w_t2->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* q_w_t3 = new TH2F("q_w_t3","Q^2 vs W for 3 Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_w_t3->GetXaxis()->SetTitle("Reconstructed Q^2 [GeV]");
    q_w_t3->GetYaxis()->SetTitle("Reconstructed W [GeV]");

    TH2F* q_w_t4 = new TH2F("q_w_t4","Q^2 vs W for 4+ Pions",NBINS_W, MIN_W, MAX_W, NBINS_W, MIN_W, MAX_W);
    q_w_t4->GetXaxis()->SetTitle("Reconstructed Q^2 [GeV]");
    q_w_t4->GetYaxis()->SetTitle("Reconstructed W [GeV]");
