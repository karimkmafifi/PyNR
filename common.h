#ifndef COMMON_H
#define COMMON_H

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#include <QtCore/QSignalMapper>
#include <QtCore/QString>
#include <QFileInfo>
#include <QThread>
#include <QMutex>
#include <qtextstream.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <thread>
#include <ctime>
#include <stack>
#include <includelibs/glm/glm/glm.hpp>
#include <includelibs/glm/glm/gtc/matrix_transform.hpp>
#include <includelibs/glm/glm/gtc/type_ptr.hpp>
#include <math.h>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <includelibs/KHR/khrplatform.h>

struct Atom;
struct Bond;
struct Ring;

const glm::vec3 zero_vec(0, 0, 0);
const float max_float = (std::numeric_limits<float>::max)();
const unsigned max_unsigned = (std::numeric_limits<unsigned int>::max)();
const int max_int = (std::numeric_limits<int>::max)();
const float epsilon_float = std::numeric_limits<float>::epsilon();
const glm::vec3 max_vec(max_float, max_float, max_float);

const float APPROX_ZERO = 1.0E-6;

const float pi = float(3.1415926535897931);
const float PI_halved = pi / 2.;

struct Error_report
{
    std::string errorText;
    Error_report(const std::string& _error) { errorText = _error; }
};

static const std::string uint_ = "unsigned int";
static const std::string uint2_ = "unsigned";
static const std::string int_ = "int";
static const std::string double_ = "double";
static const std::string long_ = "long";
static const std::string float_ = "float";

class Common {
public:

    template<typename T>
    static inline T sqr(T x)
    {
        return x*x;
    }

    static float vec_distance_sqr(const glm::vec3& a, const glm::vec3& b) {
        return sqr(a.x - b.x) + \
            sqr(a.y - b.y) + \
            sqr(a.z - b.z);
    }

    static float dot_product(const glm::vec3& a, const glm::vec3& b) {
        return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    }

    static float length(const glm::vec3& a) {
        return std::sqrt(dot_product(a, a));
    }

     static float angle_radian(const glm::vec3& a, const glm::vec3& b) {
        return acos(dot_product(a, b) / (length(a)*length(b)));
    }

     template<typename tt>
     static inline tt string_to(std::string str)
     {
         try
         {
             tt result;
             std::stringstream ss;
             ss << str;
             ss >> result;

             return result;
         }
         catch (Error_report& err)
         {
             throw Error_report(err.errorText + "Error during casting string to another datatype");
         }
     }

    template<typename T>
    static T checked_convert_substring(const std::string& str, int i, int j, const std::string& dest_nature) {

        if (!(i >= 1)) { throw Error_report("Internal Error"); }
        if (!(i <= j + 1)) { throw Error_report("Internal Error"); };

        if (j > str.size()) { throw Error_report("The line is too short"); }

        // omit leading whitespace
        while (i <= j && isspace(str[i - 1]))
        {
            ++i;
        }

        const std::string substr = str.substr(i - 1, j - i + 1);

        try {
            return Common::string_to<T>(substr);
        }
        catch (Error_report& err) {
            throw Error_report(err.errorText + " (" + substr + ") " + "is not a valid " + dest_nature);
        }
    }

    static std::string omit_whitespace(const std::string& str, int i, int j) {

        if (i < 1) { i = 1; }
        if (j < i - 1) { j = i - 1; }// i >= 1
        if (j > str.size()) { j = str.size(); }

        // omit leading whitespace
        while (i <= j && isspace(str[i - 1]))
        {
            ++i;
        }

        // omit trailing whitespace
        while (i <= j && isspace(str[j - 1]))
        {
            --j;
        }

        if (!(i - 1 < str.size())) { throw Error_report("Internal Error"); };
        if (!(j - i + 1 < str.size())) { throw Error_report("Internal Error"); };

        return str.substr(i - 1, j - i + 1);
    }

    static void split(const std::string &s, char delim, std::vector<std::string> &elems)
    {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim))
        {
            elems.push_back(item);
        }
    }

    static bool contains(const std::string &s, const std::string &word)
    {
        int found = s.find(word);

        if (found != std::string::npos)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    template<typename T>
    static bool has(const std::vector<T>& v, const T& element) {
        return std::find(v.begin(), v.end(), element) != v.end();
    }

    static void rmvBackslash(std::string &input)
    {
        for (int i = 0; i < input.size(); ++i)
        {
            if (input[i] == '\\')
            {
                input[i] = '/';
            }
        }
    }

    static void rmvWhiteSpaces(std::string &input)
    {
        input.erase(remove_if(input.begin(), input.end(), isspace), input.end());
    }

    static inline void leftTrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(isspace))));
    }

    static inline void rightTrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(isspace))).base(), s.end());
    }

    static inline void bothTrim(std::string &s) {
        leftTrim(s);
        rightTrim(s);
    }

    static inline bool substring_is_blank(const std::string& str, int i, int j) { // indexes are 1-based, the substring should be non-null
        if (i < 1 || i > j + 1 || j > str.size()) { throw Error_report("Error checking if substring is blank"); }
        for (int k = i - 1; k < j; k++)
        {
            if (!(isspace(str[k])))
            {
                return false;
            }
        }
        return true;
    }

    static inline bool starts_with(const std::string& str, const std::string& start) {
        return str.size() >= start.size() && str.substr(0, start.size()) == start;
    }

    static inline std::string get_file_name(const std::string& str)
    {
        std::vector<std::string> elems;
        Common::split(str, '/', elems);
        return elems[elems.size() - 1];
    }

    template<typename T>
    static int vector_append(std::vector<T>& x, const std::vector<T>& y) { // return old size
        int old_size = x.size();
        x.insert(x.end(), y.begin(), y.end());
        return old_size;
    }

    static bool is_number(const std::string& s)
    {
        std::string::const_iterator it = s.begin();
        while (it != s.end() && isdigit(*it)) ++it;
        return !s.empty() && it == s.end();
    }
};

#endif
