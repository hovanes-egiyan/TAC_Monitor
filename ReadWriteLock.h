/*
 * ReadWriteLock.h
 *
 *  Created on: Apr 25, 2017
 *      Author: hovanes
 */

#ifndef READWRITELOCK_H_
#define READWRITELOCK_H_

extern "C" {
#include <pthread.h>
}
;

class ReadLock {
protected:
	pthread_rwlock_t& m_;
public:
	//! Construct and lock mutex
	ReadLock(pthread_rwlock_t& m) :
			m_(m) {
		pthread_rwlock_rdlock(&m_);
	}
	//! Destruct and unlock mutex
	virtual ~ReadLock() {
		pthread_rwlock_unlock(&m_);
	}
	//! Explicitly lock mutex, need to make sure you do not get into racing condition
	virtual int lock() {
		return pthread_rwlock_rdlock(&m_);
	}
	//! Explicitly unlock existing mutex
	virtual int unlock() {
		return pthread_rwlock_unlock(&m_);
	}
	pthread_rwlock_t& getLock() {
		return m_;
	}
};

class WriteLock {
protected:
	pthread_rwlock_t& m_;
public:
	//! Construct and lock mutex
	WriteLock(pthread_rwlock_t& m) :
			m_(m) {
		pthread_rwlock_wrlock(&m_);
	}
	//! Destruct and unlock mutex
	virtual ~WriteLock() {
		pthread_rwlock_unlock(&m_);
	}
	//! Explicitly lock mutex, need to make sure you do not get into racing condition
	virtual int lock() {
		return pthread_rwlock_wrlock(&m_);
	}
	//! Explicitly unlock existing mutex
	virtual int unlock() {
		return pthread_rwlock_unlock(&m_);
	}
	pthread_rwlock_t& getLock() {
		return m_;
	}
};

#endif /* READWRITELOCK_H_ */
