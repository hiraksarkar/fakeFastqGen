#ifndef __FILESYSTEM_OP_HPP__
#define __FILESYSTEM_OP_HPP__

    // Taken from
    // http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
    bool FileExists(const char* path);
    // Taken from
    // http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
    bool DirExists(const char* path);
    void MakeDir(const char* path);

#endif //__FILESYSTEM_OP_HPP__
