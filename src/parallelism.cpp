/* mlev-hiv-model - A program for running multi-level HIV-1 simulations
 * Copyright (C) 2017 Christiaan H. van Dorp <chvandorp@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "parallelism.hpp"


std::mutex ThreadSafeCout::coutMutex;
// the cout_mutex is shared between all ThreadSafeCout objects
ThreadSafeCout::ThreadSafeCout() {	// upon construction, lock the coutMutex
	coutMutex.lock();
}
ThreadSafeCout::~ThreadSafeCout() {	// upon destruction, unlock the coutMutex
	coutMutex.unlock();
}


// methods for Job

Job::Job() {
	concurrent = false;
	finished = false;
	priority = 0;
}
Job::~Job() { /* empty */ }

bool Job::isFinished() const {
	std::lock_guard<std::mutex> lock(finishedMutex);
	return finished;
}
bool Job::waitForJobToFinish() const {
	if ( !concurrent ) {
		return false;
	}
	// make sure Job is handled quickly
	prioritize();
	// wait for the condition variable while finished is false
	std::unique_lock<std::mutex> lock(finishedMutex);
	finishedConditionalVar.wait(lock, [this]{ return finished; });
	return true;
}
void Job::setFinished() {
	std::lock_guard<std::mutex> lock(finishedMutex);
	finished = true;
	finishedConditionalVar.notify_all();
}

void Job::prioritize() const {
	/* increase the priority of the Job
	 * only has an effect when the job is in the jobs list
	 */
	std::lock_guard<std::mutex> lock(priorityMutex);
	priority++;
}

/* aux function for sorting the jobs list (friend of Job)
 * NB: this function is NOT THREAD SAFE. lock jobsMutex before sorting
 */
bool compareJobPtrs(const Job* leftJob, const Job* rightJob) {
	// a Job with a higher priority comes first
	return leftJob->priority > rightJob->priority;
}

// methods for WorkerPool

WorkerPool::WorkerPool() {
	noMoreJobsFlag = true;
	jobsyncCounter = 0;
}
WorkerPool::WorkerPool(int numOfWorkerThreads, Rng & rng) {
	initWorkerPool(numOfWorkerThreads, rng);
}

WorkerPool::~WorkerPool() {
	// empty TODO: should it be?
}

// auxiliary function used to start thread with a class method
void WorkerPool::workerThreadEntryFunc(WorkerArg warg) {
	warg.wp->workerThreadEntry(warg.seed);
}

bool WorkerPool::initWorkerPool(int numOfWorkerThreads, Rng & rng) {
	noMoreJobsFlag = false;
	jobsyncCounter = 0;
	if ( numOfWorkerThreads < 1 ) return false;
	// else...
  for ( int i = 0; i < numOfWorkerThreads; ++i ) {
		WorkerArg wa(this, rng()); // rng() produces a new seed
		workerThreads.push_back(std::thread(workerThreadEntryFunc, wa));
		// thread is moved to list, wa is passed by value
	}
	return true;
}

void WorkerPool::waitForWorkerPoolToExit() {
	sendNoMoreJobsSignal(); // otherwise the workers will not return
	for ( auto & thread : workerThreads ) thread.join();
}

void WorkerPool::syncWorkerThreads(bool verbose) {
	std::unique_lock<std::mutex> lock(jobsyncMutex);
	EtaEstimator eta(jobsyncCounter);
	while ( jobsyncCounter > 0 ) {
		jobsyncConditionalVar.wait(lock);
		eta.update();
		if ( verbose ) std::cout << "\033[2K\rETA: " << eta << std::flush;
	}
}

void WorkerPool::sendNoMoreJobsSignal() {
	std::lock_guard<std::mutex> lock(jobsMutex);
	noMoreJobsFlag = true;
	jobsConditionalVar.notify_all(); // wake all waiting threads
}

void WorkerPool::addNewJob(Job* job) {
	if ( job == nullptr ) return; // don't add NULL jobs: can cause a stall

	/* HACK: only when the job is added to the workerpool,
	 * waitForJobToFinish should halt. WARNING: Job::concurrent should only be
	 * handled by the main thread.
	 */
	job->concurrent = true;

	jobsyncMutex.lock();
	jobsyncCounter++;
	jobsyncMutex.unlock();

	jobsMutex.lock();
	jobs.push_back(job);
	jobsConditionalVar.notify_one(); // wake a waiting thread
	jobsMutex.unlock();
}

void WorkerPool::workerThreadEntry(unsigned long seed) {
	// init a thread-local RNG
	Rng rng(seed);
	// start the routine
	Job* myJob = nullptr;
	bool lookForAJob = true;
	while ( lookForAJob ) {
		if ( myJob != nullptr ) {
			myJob->execute(rng); // execute the job, true means concurrent execution
			myJob->setFinished();
			myJob = nullptr; // reset my job
			jobsyncMutex.lock();
			jobsyncCounter--;
			// send signal: syncWorkerThreads() could be listening...
			jobsyncConditionalVar.notify_one();
			jobsyncMutex.unlock();
		}
		else { // try to pop a Job from jobs
			std::unique_lock<std::mutex> jobsLock(jobsMutex);
			if ( !jobs.empty() ) {
				/* sort the jobs to get the job with the highest priority.
				 * The resulting order of equivalent elements is stable: i.e.,
				 * equivalent elements preserve the relative order that
				 * they had before the call.
				 */
				jobs.sort(compareJobPtrs);
				myJob = jobs.front();
				jobs.pop_front();
			} else if ( noMoreJobsFlag ) {
				// noMoreJobsFlag is guarded by the jobsMutex
				lookForAJob = false;
			} else {
				// wait for a new Job to be pushed to jobs
				jobsConditionalVar.wait(jobsLock);
			} // if, else if, else
		} // jobsLock is destroyed...
	} // while ( lookForAJob )
}
