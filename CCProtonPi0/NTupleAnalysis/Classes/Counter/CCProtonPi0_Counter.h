#ifndef CCProtonPi0_Counter_h
#define CCProtonPi0_Counter_h

#include<iostream>
#include<string>

class CCProtonPi0_Counter 
{
    public:
        CCProtonPi0_Counter();
        void setName(std::string inputName);
        void increment(); 
        void print();
        double getCount();
        std::string getName();

        bool isCounted;
    
    private:
        std::string name;
        double count;
};

#endif

