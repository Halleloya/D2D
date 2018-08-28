#include <iostream>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <direct.h>
using namespace std;

int main()
{
//生成参数设定
int deviceNum = 20;  //设备数
int sampleNum = 10;  //生成样本数目

///////////////////////////////////////////////////
    srand(time(NULL));
    ofstream fout;
    char str[500];
    int fz;

    for (int i=1; i<=sampleNum ;i++)
    {
        sprintf(str,"sample%d",deviceNum);     //创建文件夹
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
