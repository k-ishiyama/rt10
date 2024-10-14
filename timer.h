#pragma once
#include <chrono>

struct Timer
{
	long long GetElapsedTime() const
	{
		return std::chrono::duration_cast<std::chrono::microseconds>(GetTimeSinceEpoch() - mStartTime).count();
	}

	float GetElapsedTimeSecond() const
	{
		return 1e-6f * GetElapsedTime();
	}

	Timer()
	{
		mStartTime = GetTimeSinceEpoch();
	}

	~Timer() = default;
private:
	inline std::chrono::microseconds GetTimeSinceEpoch() const
	{
		return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now().time_since_epoch());
	}
	std::chrono::microseconds mStartTime;	// [us]
};