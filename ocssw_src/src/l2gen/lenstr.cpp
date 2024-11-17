#include <string>
#include <stdint.h>
int32_t lenstr(const std::string & str)
{
      size_t found = str.find(' ');
      if (found != std::string::npos)
            return static_cast<int32_t>(found);
      else
            return static_cast<int32_t>(str.size());
}
extern "C" int32_t lenstr_(char * str)
{
      return lenstr(str);
}