#include "timing.h"

using namespace std::chrono;

//----------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------
timing::timing(std::vector<std::string> names)
{
    m_algorithms = names;
	length = names.size();
	currentStep.reserve(length);
}


//----------------------------------------------------------------------
// Deconstructor
//----------------------------------------------------------------------
timing::~timing()
{
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

void timing::start() {
	t1 = high_resolution_clock::now();
}


void timing::stop() {
	t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	currentStep.push_back(time_span.count());
}


void timing::endStep() {
	m_timing.push_back(currentStep);
	currentStep.clear();
}


void timing::writeTiming(std::ofstream *timingOut) {
	if (timingOut->is_open()) {
		for (int i = 0; i < length; i++) {
			double tim_j = 0;
			for (int j = 0; j < m_timing.size(); j++) {
				tim_j += m_timing[j][i];
			}
			tim_j = tim_j / m_timing.size();
			double sig_j = 0;
			for (int j = 0; j < m_timing.size(); j++) {
				sig_j += (tim_j - m_timing[j][i]) * (tim_j - m_timing[j][i]);
			}
			sig_j = sig_j / m_timing.size();
			sig_j = sqrt(sig_j);

			(*timingOut) << m_algorithms[i] << ", " << tim_j << ", " << sig_j << "\n";
		}		
	}
}
