#ifndef Par_h_seen
#define Par_h_seen

class Par {
public:
    static Par& Get();
    
    double k_ecal;
    double k_hcal;
    double k_trkr;
    double coangle;
    double clength;
    double xminevis;
    double uminevis;
    double vminevis;

    void Print();

private:
    Par() {}
    Par(const Par&) {}
    Par& operator=(const Par&) { return *this;}
    ~Par() {}
};

Par& Constants();

#endif
