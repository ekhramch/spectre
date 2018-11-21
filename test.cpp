#include <Magick++.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <windows.h>
#include <tchar.h> 
#include <stdio.h>
#include <strsafe.h>
#pragma comment(lib, "User32.lib")
#include <vector>
#include "cub_spline.h"
#include <math.h>
#include <conio.h>

using namespace std;

using namespace Magick;

struct E
{
	int groups_count;

	int group_number;

	vector<int> et_in_group;

	E(int n, int m):et_in_group(n,0)
	{
		int tmp=0;

		groups_count=n;

		group_number=m;

		if(n>1)
		{
			for(int i=0; i<n; i++)
			{
				cout << endl << "Введите количество эталонов в группе №" << i+1 << endl;

				cin.clear();
				cin.sync();

				cin >> tmp;

				if(cin.fail())
				{
					cout << "Введено некорректное значение!" << endl;
					i--;				
				}

				et_in_group[i]=tmp;
			}
		}
	}

	E():et_in_group(1,0)
	{
		groups_count=1;
		group_number=0;
	}
};

int img_med(vector<double> &_med, vector<string> &_file_names, Image &_image, string folder, string flags);

int evaluation(vector<double> &_probes, vector<string> &_file_names_probes, vector<double> &_conc_probes, Image &_image, string folder, E &_etalon, string flags);

int interpolation_spline(vector<double> &_etals, vector<double> &_conc_etals, vector<double> &_probes, vector<double> &_conc_probes);

int interpolation_line(vector<double> &_etals, vector<double> &_conc_etals, vector<double> &_probes, vector<double> &_conc_probes);

int interpolation(vector<double> &_etals, vector<double> &_conc_etals, vector<double> &_probes, vector<double> &_conc_probes, string flag);

int file_control(vector<double> &_probes, vector<string> &_file_names_probes, vector<double> &_conc_probes, string filename);

string set_flag();

int parsing();

int main( int argc, char ** argv)
{
	vector<string> file_names_etals;
	vector<double> etals;
	vector<double> conc_etals;
	vector<double> med_etals;

	vector<string> file_names_probes;
	vector<string> file_names_probes_med;
	vector<double> probes;
	vector<double> conc_probes;
	vector<double> med_probes;
	
	setlocale(0, "rus");
	
	InitializeMagick(*argv);
	Image image;

	string flags=set_flag();

	string int_method = "ols";

	int count, et_num;
	
	cout << endl << "Введите количество эталонных групп" << endl;

	do
	{
		cin.clear();
		cin.sync();
		cin >> count;
		if(cin.fail())
			cout << "Введите корректное значение!" << endl;
	}
	while(cin.fail());

	if(count>1)
	{
		cout << endl << "Введите номер целевой эталонной группы" << endl;

		do
		{
			cin.clear();
			cin.sync();
			cin >> et_num;
			if(cin.fail())
				cout << "Введите корректное значение!" << endl;
		}
		while(cin.fail());

		et_num--;
	}
	else
		et_num=0;

	E etalon(count,et_num);

	evaluation(etals, file_names_etals, conc_etals, image, "etals", etalon, flags);

	evaluation(probes, file_names_probes, conc_probes, image, "probes", etalon, flags);

	interpolation(etals, conc_etals, probes, conc_probes, int_method);

	file_control(probes, file_names_probes, conc_probes, "control-probes.txt");

	parsing();

	return 0;
}


int interpolation_spline(vector<double> &_etals, vector<double> &_conc_etals, vector<double> &_probes, vector<double> &_conc_probes)
{
	cubic_spline S;

	for(int i=0; i<_conc_etals.size();i++)
		_conc_etals[i] = log10(_conc_etals[i]);

	S.build_spline(_etals.data(), _conc_etals.data(), _conc_etals.size());

	for(int i=0; i<_probes.size(); i++)
	{
		_conc_probes[i] = S.f(_probes[i]);
		_conc_probes[i] = pow(10.0,_conc_probes[i]);
	}

	ofstream odata("data_fin_tmp.txt");
	ostream_iterator<double> iout(odata, "\n");
	copy(_conc_probes.begin(), _conc_probes.end(), iout);
	odata.close();

	return 0;
}

int interpolation_line(vector<double> &_etals, vector<double> &_conc_etals, vector<double> &_probes, vector<double> &_conc_probes)
{
	double S1=0.0, S2=0.0, S3=0.0, S4=0.0;

	for(int i=0; i<_conc_etals.size();i++)
		_conc_etals[i] = log10(_conc_etals[i]);

	if(_etals.size()==_conc_etals.size())
	{
		for(int i=0; i<_etals.size(); i++)
		{
			S1 = S1 + _etals[i];
			S3 = S3 + _etals[i]*_etals[i];
		}

		for(int i=0; i<_conc_etals.size(); i++)
		{
			S2 = S2 + _conc_etals[i];
			S4 = S4 + _conc_etals[i]*_etals[i];
		}
	}
	else
	{
		cout << "ERROR! NOT ENOUGH ETALON CONCENTRATION! HALTED!" << endl;
		exit(-1);
	}

	double A=0.0, B=0.0, N=static_cast<double>(_conc_etals.size());

	A = (N*S4 - S1*S2)/(N*S3 - S1*S1);
	B = (S2 - S1*A)/N;

	if(_probes.size()==_conc_probes.size())
	{
		for(int i=0; i<_probes.size(); i++)
			_conc_probes[i] = A*_probes[i] + B;
	}
	else
	{
		cout << "ERROR!" << endl;
		exit(-1);
	}

	for(int i=0; i<_probes.size(); i++)
		_conc_probes[i] = pow(10.0,_conc_probes[i]);

	ofstream odata("data_fin_tmp.txt");
	ostream_iterator<double> iout(odata, "\n");
	copy(_conc_probes.begin(), _conc_probes.end(), iout);
	odata.close();

	return 0;
}

int file_control(vector<double> &_probes, vector<string> &_file_names_probes, vector<double> &_conc_probes, string filename)
{
	FILE *data;
	data=fopen(filename.c_str(), "w" );
	for(int i=0; i<_probes.size(); i++)
		fprintf(data,"%s %f %f\n", _file_names_probes[i].c_str(), _probes[i], _conc_probes[i]);
	fclose(data);

	return 0;
}

int img_med(vector<double> &_med, vector<string> &_file_names, Image &_image, string folder, string flags)
{
	WIN32_FIND_DATA ffd;	
	TCHAR szDir[MAX_PATH];
	HANDLE hFind = INVALID_HANDLE_VALUE;

	string tmp = "./" + folder;
	StringCchCopy(szDir, MAX_PATH, tmp.c_str());
	StringCchCat(szDir, MAX_PATH, TEXT("\\*"));
	hFind = FindFirstFile(szDir, &ffd);

	PixelPacket *pixels;

	do
	{
		if (!(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
		{
			string s(ffd.cFileName);
			_file_names.push_back(s);
			_med.push_back(0.0);
		}
	}
	while (FindNextFile(hFind, &ffd) != 0);

	for(int iter=0; iter<_file_names.size(); iter++)
	{
		_image.read(".\\"+folder+"\\"+_file_names[iter]);

		int w = _image.baseColumns();
		int h = _image.baseRows();
		pixels = _image.getPixels(0, 0, w, h);
		
		for(int i=0; i<w*h; i++)
			_med[iter]+= pixels[i].blue*0.114 + pixels[i].green* 0.587 + pixels[i].red* 0.299;

		_med[iter] = _med[iter]/( static_cast<double>(h) * static_cast<double>(w) );
	}

	return 0;
}

int evaluation(vector<double> &arr, vector<string> &_file_names, vector<double> &conc_arr, Image &_image, string folder, E &_etalon, string flags="iml")
{
	WIN32_FIND_DATA ffd;	
	TCHAR szDir[MAX_PATH];
	HANDLE hFind = INVALID_HANDLE_VALUE;

	string tmp = "./" + folder;
	StringCchCopy(szDir, MAX_PATH, tmp.c_str());
	StringCchCat(szDir, MAX_PATH, TEXT("\\*"));
	hFind = FindFirstFile(szDir, &ffd);

	PixelPacket *pixels;


	/*medians of images evaluation*/
	vector<string> file_names_med_tmp;
	vector<double> med_tmp;

	vector<double> med;
	vector<string> file_names_med;


	if(folder=="etals")
	{
		img_med(med_tmp, file_names_med_tmp, _image, "full_etals", flags);

		if(_etalon.groups_count==1)
		{
			file_names_med=file_names_med_tmp;
			
			med=med_tmp;
		}
		else
		{
			int j=0;
			
			for(int i=0; i<_etalon.group_number; i++)
				j+=_etalon.et_in_group[i];

			file_names_med.insert(file_names_med.begin(), file_names_med_tmp.begin() + j, file_names_med_tmp.begin() + j + _etalon.et_in_group[_etalon.group_number]);

			med.insert(med.begin(), med_tmp.begin() + j, med_tmp.begin() + j + _etalon.et_in_group[_etalon.group_number]);
		}
	}
	else
	{
		img_med(med_tmp, file_names_med_tmp, _image, "full_probes", flags);

		file_names_med=file_names_med_tmp;

		med=med_tmp;
	}

	/*end of medians of images evaluation*/

	/*concentration block*/
	do
	{
		if (!(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
		{
			string s(ffd.cFileName);
			_file_names.push_back(s);
			arr.push_back(0.0);
			conc_arr.push_back(0.0);
		}
	}
	while (FindNextFile(hFind, &ffd) != 0);

	if(_file_names.size()!=file_names_med.size())
	{
		cout << "Количество полных изображений должно совпадать с количеством изображений с линией! Нажмите любую клавишу." << endl;
		_getch();
		exit(-1);
	}

	if(folder=="etals")
	{
		ifstream fdata("data_conc.txt");
		istream_iterator<double> idata(fdata);
		istream_iterator<double> eof;
		if(idata==eof)
		{
			cout << "Файл концентраций data_conc.txt не содержит значений. Нажмите любую клавишу." << endl;
			_getch();
			exit(-1);
		}
		vector<double> tmp_arr(idata, eof);
		if(tmp_arr.size()!=conc_arr.size())
		{
			cout << "Недостаточное количество файлов эталонов. Количество эталонных изображений должно совпадать с количеством значений в data_conc.txt! Нажмите любую клавишу." 
				<< endl;
			_getch();
			exit(-1);
		}
		else
			conc_arr = tmp_arr;

		for(int i=0; i<conc_arr.size(); i++)
		{
			if(conc_arr[i]==0)
			{
				cout << "Недостаточное количество эталонов в data_conc.txt! Нажмите любую клавишу." 
					<< endl;
				_getch();
				exit(-1);
			}
		}
		
		if(_etalon.groups_count>1)
		{
			int j=0;

			for(int i=0; i<_etalon.group_number; i++)
				j+=_etalon.et_in_group[i];

			vector<string> names_tmp(_file_names.begin() + j, _file_names.begin() + j + _etalon.et_in_group[_etalon.group_number]);
			vector<double> arr_tmp(arr.begin() + j, arr.begin() + j + _etalon.et_in_group[_etalon.group_number]);
			vector<double> conc_tmp(conc_arr.begin() + j, conc_arr.begin() + j + _etalon.et_in_group[_etalon.group_number]);

			_file_names = names_tmp;
			arr = arr_tmp;
			conc_arr = conc_tmp;
		}
	}

	for(int iter=0; iter<_file_names.size(); iter++)
	{
		_image.read(".\\"+folder+"\\"+_file_names[iter]);

		int w = _image.baseColumns();
		int h = _image.baseRows();
		pixels = _image.getPixels(0, 0, w, h);

		for(int i=0; i<w*h; i++)
			arr[iter]+= pixels[i].blue*0.114 + pixels[i].green* 0.587 + pixels[i].red* 0.299;


		arr[iter] = arr[iter]/( static_cast<double>(h) * static_cast<double>(w) );

		if(flags[1]=='w')
			arr[iter] = 1.0 - (arr[iter] - med[iter])/65536.0;
		else
		{
			if(flags[1]=='m')
				arr[iter] = arr[iter]/med[iter];
			else
				arr[iter] = 1.0 - arr[iter]/65536.0;
		}

		/*if(flags[2]=='l')
			if(arr[iter]>0)
				arr[iter] = log10(arr[iter]);
			else
			{
				cout << "Логарифмическая шкала не может быть установлена для неположительных величин. Попробуйте другую шкалу или обратитесь к разработчику." << endl;
				_getch();
			}
		else
			if(flags[2]=='p')
				arr[iter] = pow(10,arr[iter]);*/

		if(folder=="etals")
			file_control(arr, _file_names, conc_arr, "control-etals.txt");
	}

	return 0;
}

int interpolation(vector<double> &_etals, vector<double> &_conc_etals, vector<double> &_probes, vector<double> &_conc_probes, string flag="ols")
{
	if(flag=="spline")
		interpolation_spline(_etals, _conc_etals, _probes, _conc_probes);
	else
		interpolation_line(_etals, _conc_etals, _probes, _conc_probes);

	return 0;
}

int parsing()
{
	vector<string> pars;

	ifstream fdata("data_fin_tmp.txt");
	istream_iterator<string> idata(fdata);
	istream_iterator<string> eof;
	pars.insert(pars.begin(),idata, eof);
	fdata.close();


	for(int i=0; i<pars.size();i++)
		for(int j=0; j<pars[i].size();j++)
				if(pars[i][j]=='.')
					pars[i][j]=',';

	ofstream odata("data_fin.txt");
	ostream_iterator<string> iout(odata, "\n");
	copy(pars.begin(), pars.end(), iout);

	remove("data_fin_tmp.txt");

	return 0;
}

string set_flag()
{
	string flag="i";

	cout << "Выберите метод осреднения: w - чистым белым с вычитанием фона" << endl << "			   m - средним фоном рисунка без вычитания фона" 
		<< endl << "			   s - чистым белым без вычитания фона" << endl ;

	char c;
	cin >> c;

	if(c=='w')
	{
		flag+= "w";
		cout << "Выбран вариант: чистым белым с вычитанием фона" << endl << endl;
	}
	else
	{
		if(c=='s')
		{
			flag+= "s";
			cout << "Выбран вариант: чистым белым без вычитания фона" << endl << endl;
		}
		else
		{
			flag+= "m";
			cout << "Выбран вариант: средним фоном рисунка без вычитания фона" << endl << endl;
		}
	}

	flag+= "l";
	
	/*cout << "Выберите шкалу почернения: n - стандартная" << endl << "			   l - логарифмическая" 
		<< endl << "			   p - степень 10" << endl;
	
	cin >> c;

	if(c=='l')
	{
		flag+= "l";
		cout << "Логарифмическая шкала" << endl << endl;
	}
	else
	{
		if(c=='p')
		{
			flag+= "p";
			cout << "Степенная шкала" << endl << endl;
		}
		else
		{
			flag+= "n";
			cout << "Стандартная шкала" << endl << endl;
		}
	}*/

	return flag;
}