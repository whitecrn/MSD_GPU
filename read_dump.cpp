#include "data.h"

void read_dump
(vector<atom_old>& R_old,int &atom_number,int &step,double &a,double &b,double &c)
{    
    double x_b[6];
    int total_lines=0;
    string line,word;
    int lines;
    ifstream file;
    file.open(file_name);
    if (!file.is_open())
    {
        cout << "Can't open dump.atom file." << endl;
        exit(1); 
    }

//comfirm atom_number:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (int i=0;i<3;i++)
    {
        getline(file,line);
    }
    file >> atom_number;   //read the total atom number
    getline(file,line);
    getline(file,line);
    for (int i=0;i<3;i++)
    {
        getline(file,line);
        istringstream WORDS(line);
        for (int j=0;j<2;j++)
        {
            if (WORDS >> word)
            {
                if (j==0)
            {
                x_b[i*2+j]=stod(word);
            }
            else if (j==1)
            {
                x_b[i*2+j]=stod(word);
            }
            }
        }
    }
    a=x_b[1]-x_b[0];
    b=x_b[3]-x_b[2];
    c=x_b[5]-x_b[4];
    getline(file,line);
    istringstream WORDS(line);
    if (WORDS >> word)
    {
        if (word=="xy")
            lines=10;
        else
            lines=9;
    }
    file.clear();
    file.seekg(0,ios::beg);
//comfirm steps:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    while(getline(file,line))
    {
        ++total_lines;
    }
    step=total_lines/(lines+atom_number);
    file.clear();
    file.seekg(0,ios::beg);
//read xyz:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (int i=0;i<step;i++)
    {
        for (int j=0;j<lines;j++)
        {
            getline(file,line);
        }

        for (int j=0;j<atom_number;j++)
        {
            getline(file,line);
            istringstream words(line);
            int indices,type;
            double x,y,z;
            for (int k=0;k<5;k++)
            {
                if (words >> word)
                {
                    if (k==0)
                    {
                        indices=stoi(word);
                    }
                    else if (k==1)
                    {
                        type=stoi(word)-1;
                    }
                    else if (k==2)
                    {
                        x=stod(word);
                    }
                    else if (k==3)
                    {
                        y=stod(word);
                    }
                    else if (k==4)
                    {
                        z=stod(word);
                    }
                }
            }
            R_old.push_back(atom_old(i,indices,type,x,y,z));
        }
        if (i%static_cast<int>(round(0.2*step))==0)
        {
            cout << i << " is read. " << endl;
        }
    }
    sort(R_old.begin(),R_old.end());
}