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

#ifndef PARALLELISM_HPP
#define PARALLELISM_HPP

#include <thread>
#include <mutex>
#include <condition_variable>
#include <list>
#include <iostream>

#include "eta.hpp"
#include "rng.hpp"
// todo: think of a better way to have multiple Rng's

struct ThreadSafeCout {
	/** ThreadSafeCout is supposed to act as std::cout,
	 * but it lets a thread print un-interrupted lines by locking a
	 * static mutex. ThreadSafeCout objects are supposed to remain
	 * unnamed. The mutex is locked upon creation and unlocked upon
	 * destruction. e.g.
	 * ThreadSafeCout() << "hello from thread " << pthread_self() << "\n";
	 * locks the shared mutex, prints a line, and unlocks the mutex after
	 * the semicolon.
	 * TODO/FIXME: inherit basic_ostream. ThreadSafeCout can not handle
	 * input like std::endl and std::flush. Therefore std::flush is passed
	 * to std::cout after each call of operator<<
	 */
	static std::mutex coutMutex; // shared mutex
	ThreadSafeCout(); // lock the mutex
	~ThreadSafeCout(); // unlock the mutex
	template<typename T>
	ThreadSafeCout & operator<<(T x) { // pass x to std::cout
		std::cout << x << std::flush;
		return *this; // return a reference of the ThreadSafeCout object
	}
};

/* todo:
 * (1) when jobs are easy, it might be better to pass multiple
 * jobs at once (i.e. in serial) to the workers. Maybe this
 * could be determined automatically?
 * (2) add a method addJobs that accepts ranges.
 * (3) add a finished_jobs list to worker pool so that they can be properly
 * deleted.
 */

/* the interface Job is used by WorkerPool, which contains
 * a list of Jobs. Classes that inherit Job, can be handled
 * by WorkerPool, and must implement the method "execute"
 */
class Job {
	/* WorkerPool::workerThreadEntry needs to call setFinished()
	 * after Job::execute returns. But, I don't want anybody else to
	 * call setFinished! todo: better method?
	 */
	friend class WorkerPool;
	friend bool compareJobPtrs(const Job* , const Job* );
public:
	Job(); // init private members
	virtual ~Job(); // obligatory virtual destructor
	bool isFinished() const; // NB can return 'false' incorrectly.
	/* HACK: waitForJobToFinish immediately returns false whenever the Job has not
	 * been passed to WorkerPool::addJobs.
	 */
	bool waitForJobToFinish() const;
protected:
	/* WorkerPool::workerThreadEntry calls this method.
	 * TODO Replace Rng with general type (template?)
	 */
	virtual bool execute(Rng & )=0;
private: // derived classes don't need to mess with these members...
	bool concurrent;
	bool finished;
	// add the specifier mutable so that it can be used in const methods
	mutable std::mutex finishedMutex;
	mutable std::condition_variable finishedConditionalVar;
	mutable int priority;
	mutable std::mutex priorityMutex;
	void setFinished();
	void prioritize() const;
};

class WorkerPool {
public:
	/* constuctor */
	WorkerPool();
	WorkerPool(int numOfWorkerThreads, Rng & rng); // calls initWorkerPool
	/* destructor */
	virtual ~WorkerPool();
	/* Initiate private members and set the number of workers.
	 * Returns true if the threads were successfully started,
	 * false if there was an error starting a thread
	 */
	bool initWorkerPool(int , Rng & );
	/* Will not return until the internal threads have exited. */
	void waitForWorkerPoolToExit();
	/* Wait for all Workers to finish their jobs. The Workers
	 * will then be waiting for new jobs. (So no joining...)
	 */
	void syncWorkerThreads(bool verbose=false);
	/* Send the worker pools a signal that there will be no more jobs.
	 * This will make them exit their WorkerThreadEntry function.
	 */
	void sendNoMoreJobsSignal();
	/* add a new job to the list of jobs.
	 */
	void addNewJob(Job* );
protected:
	// no protected members?
private:
	/* Upon creation of a thread, an argument can be passed. This argument
	 * needs to contain a pointer to this object (WorkerPoolWrapper),
	 * but we can add other stuff... In this case each worker gets a
	 * seed for a thread-local rng.
	 */
	struct WorkerArg {
		// inline constructors
		WorkerArg(WorkerPool* wp, unsigned long seed) : wp(wp), seed(seed) {}
		// the actual data...
		WorkerPool* wp; // a pointer to the workerpool
		unsigned long seed;
	};
	// private methods
	void workerThreadEntry(unsigned long ); // pass a seed
	static void workerThreadEntryFunc(WorkerArg );
	std::mutex jobsMutex;
	std::condition_variable jobsConditionalVar;
	std::mutex jobsyncMutex;
	std::condition_variable jobsyncConditionalVar;
	int jobsyncCounter;
	std::list<Job*> jobs;
	bool noMoreJobsFlag;
	std::list<std::thread> workerThreads;
	int numOfWorkerThreads;
};

#endif
