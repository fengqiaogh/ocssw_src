
#ifndef FOCS_DATFILE
#define FOCS_DATFILE

#include <functional>
#include <iterator>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace focs {
    /*!
        \class focs::DatFile
    
        \brief .dat file reader
    
        .dat files are key-value stores primarily used to store arrays of
        values.  Singular values are supported, but this implementation treats
        them as vectors with a single element.  On top of that, the indices
        given in the file are honored, so most variable indices start at 1.
    
        \section Format
        \include testdata/libfocs/datfile.dat

        \section Usage
        \snippet tests/datfile.cpp DatFile example
    */
    class DatFile {
        using Callback = std::function<bool(const std::tuple<std::string, int, double>&)>;
        public:
            /*!  \brief Default constructor, creates an empty object.  */
            DatFile() = default;
            /*!
                \brief Create a new object and load the given file upon creation.

                \param file Path of file to load
            */
            explicit DatFile(const std::string& file);
            /*!
                \brief Load a file, calling a user function for each line read

                \param file Path of file to load
                \param callback User function to call
            */
            DatFile(const std::string& file, Callback callback);

            /*!
                \brief Load dat file from path
            
                \param file Path of file to load
            
                \return true on failure, false on success
            */
            bool load(const std::string& file);
            /*!
                \brief Load dat file from input stream
            
                \param in stream from which to load data
            
                \return true on failure, false on success
            */
            bool load(std::istream& in);
            /*!
                \brief Load a single line
            
                \param line key-value line to load
            
                \return true on failure, false on success
            */
            bool load_line(const std::string& line);

            bool empty() { return kv_.empty(); };

            /*!
                \brief Iterator returning each variable read
            
                \return Beginning iterator
            */
            std::unordered_map<std::string, std::vector<double>>::const_iterator begin() const;
            /*!
                \brief End variable iterator 
            
                \return Ending iterator
            */
            std::unordered_map<std::string, std::vector<double>>::const_iterator end() const;

            /*!
                \brief Access a variable

                If the variable doesn't exist, it'll silently be created.

                \param k name of vector
            
                \return reference to variable's vector
            */
            std::vector<double>& operator[] (const std::string& k);
            /*!
                \brief Access a variable

                If the variable doesn't exist, it'll silently be created.

                \param k name of vector
            
                \return reference to variable's vector
            */
            std::vector<double>& operator[] (std::string&& k);

            /*!
                \brief Output dat file to stream
            
                \param os Target output stream
                \param kv DatFile object to print
            
                \return Given output stream, for chaining
            */
            friend std::ostream& operator<<(std::ostream& os, const DatFile& kv);

            /*!
                \brief Set the user function to call for each line
            
                \param callback User function to use
            */
            void callback(Callback callback){callback_ = callback;}
            /*!
                \brief Return the user function to call for each line
            
                \return The currently set user function
            */
            Callback callback(){return callback_;}

        private:
            std::unordered_map<std::string, std::vector<double>> kv_{};

            Callback callback_{};
    };
} // namespace focs

#endif // FOCS_DATFILE

