#include <iostream>

class Student
{
    public:
        string name;
        string major;
        double gpa;
        Student(string Name, string Major, double GPA){
            name = Name;
            major = Major;
            gpa = GPA;
        }
        bool hasHonors(){
            if(gpa>=3.5){
                return true;}
            return false;
        }

};


void test()
{
    Student student1("Jim","Business",2.4);
    Student student2("Pam","Art",3.6);

    cout<< student1.hasHonors()<<endl;
    cout<< student2.hasHonors()<<endl;
}