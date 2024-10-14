#pragma once
#include <iostream>
#include <fstream>
#include <mutex>
#include <string>
#include <cassert>
#include <format>  // C++20

class Logger
{
public:
  static Logger& GetInstance()
  {
    static Logger instance;
    return instance;
  }

  template <typename... Args>
  void Log(const std::string& inMessage, const Args&... inArguments)
  {
    std::lock_guard<std::mutex> lock(mMutex);
    std::string formatted = std::vformat(inMessage, std::make_format_args(inArguments...));
    if (!mFileStream.is_open())
    {
      assert(false && "Failed to write to log file.");
    }
    mFileStream << formatted << std::endl;
    mFileStream.flush();

    std::cout << formatted << std::endl;
  }

private:
  Logger()
  {
    mFileStream.open("log.txt", std::ios::out);
    assert(mFileStream.is_open() && "Failed to open log file.");
  }

  ~Logger()
  {
    if (mFileStream.is_open())
    {
      mFileStream.close();
    }
  }

  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;

  std::ofstream mFileStream;
  std::mutex mMutex;
};

template <typename... Args>
void Log(const std::string& inMessage, const Args&... inArguments)
{
  Logger::GetInstance().Log(inMessage, inArguments...);
}
