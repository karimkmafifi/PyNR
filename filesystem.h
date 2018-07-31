#ifndef FILESYSTEM_H
#define FILESYSTEM_H

#include "input.h"

class FileSystem
{
public:
	static Input confReader(const std::string& configPath);
	static std::vector<std::string> pdbReader(const std::string& pdbPath);
	static std::vector<std::string> shaderReader(const std::string& shaderPath);
	static bool file_exists(const std::string& name);
	static bool dir_exists(const std::string name);
	static std::string getPath(const std::string &configPath, const std::string &fileName);
	static void writeFile(const std::string &path, const std::string &data);
};

#endif
