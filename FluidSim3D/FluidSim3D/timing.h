#pragma once
#include<vector>
#include<string>
#include<ctime>
#include <chrono>
#include <ratio>
#include <iostream>
#include <fstream>

class timing
{
public:
	/*
	constructor takes in a vector of strings with the names of the algorithms
	*/
	timing(std::vector<std::string>);
	~timing();
	/*
	Start timer clock
	*/
	void start();
	/*
	Stop timer clock and store time difference
	*/
	void stop();
	/*
	Store current measurement series
	*/
	void endStep();
	/*
	Output average values of the measurement series for each algorithm
	Format: "algorithmname, avg-value"
	*/
	void writeTiming(std::ofstream *timingOut);
private:
	int length;
	std::vector<std::vector<double>> m_timing;
	std::vector<std::string> m_algorithms;
	std::vector<double> currentStep;
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
};

