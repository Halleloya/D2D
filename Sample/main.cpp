#include <iostream>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <direct.h>
using namespace std;

int main()
{
//���ɲ����趨
int deviceNum = 20;  //�豸��
int sampleNum = 10;  //����������Ŀ

///////////////////////////////////////////////////
    srand(time(NULL));
    ofstream fout;
    char str[500];
    int fz;

    for (int i=1; i<=sampleNum ;i++)
    {
        sprintf(str,"sample%d",deviceNum);     //�����ļ���
        mkdir(str);
        sprintf(str,"sample%d\\sample%d_%d.in",deviceNum,deviceNum,i);
        fout.open(str);

        fout<<deviceNum<<'\t'<<-1<<'\t'<<-1<<endl;

        for (int j=0; j<deviceNum; j++)
        {
            fz=0;
            fout<<fz;
            for (int k=0; k<deviceNum; k++)
            {
                fout <<'\t'<< (fz+=(1+rand()%9));
            }
            fout<<endl;

        }

        fout.close();
    }
    return 0;
}
