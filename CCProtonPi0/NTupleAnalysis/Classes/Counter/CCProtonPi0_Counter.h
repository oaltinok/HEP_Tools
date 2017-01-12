#ifndef CCProtonPi0_Counter_h
#define CCProtonPi0_Counter_h

#include<iostream>
#include<string>

class CCProtonPi0_Counter 
{
    public:
        CCProtonPi0_Counter();
        void setName(std::string inputName);
        void increment(double wgt = 1.0); 
        void print();
        double getCount();
        double calcPercent(double reference);
        std::string getName();

        bool isCounted;
    
    private:
        std::string name;
        double count;
};

#endif

