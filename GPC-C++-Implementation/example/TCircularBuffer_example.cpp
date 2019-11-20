
#include <iostream>
#include <TMath.h>
#include <fstream>

void test_stream() {
	using namespace std;
	TCircularBuffer buffer(10);
	for (int i = 10; i < 34; ++i) {
		buffer.push(i);
	}

	{
		std::cout << "Before writing buffer: " << std::endl;
		std::cout << buffer << std::endl;
		ofstream fout("data.txt");
		fout << buffer;
		fout.close();
		std::cout << "After loading buffer: " << std::endl;
		ifstream fin("data.txt");
		fin >> buffer;
		std::cout << buffer << std::endl;
	}{
		std::cout << "Before writing buffer: " << std::endl;
		std::cout << buffer << std::endl;
		ofstream fout("data.txt", std::ios::binary);
		fout << buffer;
		fout.close();
		std::cout << "After loading buffer: " << std::endl;
		ifstream fin("data.txt", std::ios::binary);
		fin >> buffer;
		std::cout << buffer << std::endl;
	}

}

int main(void){
	TCircularBuffer buffer(10);
	for (int i = 10; i < 30; ++i) {
		buffer.push(i);
		std::cout << buffer << std::endl << buffer.differ()<< std::endl;
	}

	{
		std::cout << "toReverse: " << std::endl;
		auto rev = buffer.toReverse();
		std::cout << "buffer = " << buffer << std::endl;
		std::cout << "rev = " << rev << std::endl;
	}

	{
		buffer.push(100);
		std::cout << "reverse: " << std::endl;
		std::cout << "buffer = " << buffer << std::endl;
		auto rev = buffer.reverse();
		std::cout << "rev = " << rev << std::endl;
	}

	TCircularBuffer buffer1(10), buffer2(10);
	for (int i = 10; i < 20; ++i) {
		buffer1.push(i);
		buffer2.push(100L - i);
	}
	buffer2.push(200L);
	std::cout << buffer1 << std::endl << buffer2 << std::endl;
	buffer1 += buffer2;
	std::cout << buffer1 << std::endl;


	std::cout << std::endl<< std::endl<< std::endl<< "Begin the addition test" << std::endl;

	{
		double datA[]{ 0.2100,0.1470,0.3430,0};
		TCircularBuffer bufferA(datA, 4);
		double datB[]{0.0900,0.0630,0.0441,0.1029};
		TCircularBuffer bufferB(datB, 4);
		std::cout << bufferA << std::endl 
			<< bufferB << std::endl 
			<< (bufferA + bufferB) << std::endl;
		std::cout << "----------------------------------------" << std::endl;
	}
	{
		double datA[]{0.2100,0.3871,0.1029,0};
		TCircularBuffer bufferA(datA, 4);
		double datB[]{0.0900,0.0630,0.0441,0.1029};
		TCircularBuffer bufferB(datB, 4);
		std::cout << bufferA << std::endl 
			<< bufferB << std::endl 
			<< (bufferA + bufferB) << std::endl;
		std::cout << "----------------------------------------" << std::endl;
	}

	{
		double datA[]{0.4501,0.1470,0.1029,0};
		TCircularBuffer bufferA(datA, 4);
		double datB[]{0.0900,0.0630,0.0441,0.1029};
		TCircularBuffer bufferB(datB, 4);
		std::cout << bufferA << std::endl 
			<< bufferB << std::endl 
			<< (bufferA + bufferB) << std::endl;
		std::cout << "----------------------------------------" << std::endl;
	}

	{
		double datA[]{0.2100,0.1470,0.1029,0};
		TCircularBuffer bufferA(datA, 4);
		double datB[]{0.1620,0.1134,0.0794,0.1853};
		TCircularBuffer bufferB(datB, 4);
		std::cout << bufferA << std::endl 
			<< bufferB << std::endl 
			<< (bufferA + bufferB) << std::endl;
		std::cout << "----------------------------------------" << std::endl;
	}

	try{
	test_stream();
	}catch(...){
		std::cerr << "Something bad happened during the I/O process" << std::endl;
		return -1;
	}
	return 0;
}
