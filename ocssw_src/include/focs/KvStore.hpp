
#ifndef FOCS_KVSTORE
#define FOCS_KVSTORE

#include "focs/FilterIterator.hpp"

#include <functional>
#include <iterator>
#include <string>
#include <unordered_map>
#include <vector>

namespace focs {
    /*!
        \class focs::KvStore

        \brief Key-value file reader (commonly called par files)

        \section Format

        A simple file with various quoting and no commands (test1.par):
        \include testdata/libfocs/kv_store/test1.par

        A key-value file with groups and another key-value file loaded first (test2.par):
        \include testdata/libfocs/kv_store/test2.par

        \section Examples

        Creating a simple map out of key-value files:
        \snippet tests/kvstore.cpp focs::KvStore example1

        Walking though a key-value file as it is being read:
        \snippet tests/kvstore.cpp focs::KvStore example2
    */
    class KvStore {
        using read_callback = std::function<bool(const std::pair<std::string, std::string>&)>;
        using group_filtered_iterator = focs::FilterIterator<std::unordered_map<std::string, std::string>::const_iterator, std::function<bool(const std::pair<std::string, std::string>&)>>;
        using group_iterator = std::vector<std::string>::const_iterator;
        using iterator = std::unordered_map<std::string, std::string>::iterator;
        using const_iterator = std::unordered_map<std::string, std::string>::const_iterator;

        public:
            /*!  \brief Default constructor, use `load(std::string file)` to read files */
            KvStore() = default;
            /***
             * @brief
             * @param kv_store supplied key-value
             */
             explicit KvStore(std::unordered_map<std::string,std::string> & kv_store);
            /*!  \brief Construct a key-value map, loading the given file upon creation

                 \param file Path to key-value store to read

                 \throw not::sure if any lines read are invalid
             */
            explicit KvStore(const std::string& file);

            /*!  \brief Construct a key-value map, with default file and processing callback

                 \param file Path to key-value store to read
                 \param callback Function to call for each value read

                 \throw not::sure if any lines read are invalid
             */
            KvStore(const std::string& file, read_callback callback);

            /*!  \brief Read key-value file, appending values to current map

                 \param file Path to key-value store to read

                 \throw not::sure if any lines read are invalid

                 \return true on failure, false on success
            */
            bool load(const std::string& file);

            /*!  \brief Read key-values from input stream, appending values to current map

                 \param in Input stream containing lines of key-value definitions

                 \throw not::sure if any lines read are invalid

                 \return true on failure, false on success
            */
            bool load(std::istream& in);

            /*!  \brief Process key-value line, appending value to current map

                 \param line Single key-value line to process

                 \throw not::sure if any lines read are invalid

                 \return true on failure, false on success
            */
            bool load_line(const std::string& line);

            /*!  \brief std::vector iterator to recurse group names

                 \return Start iterator
            */
            group_iterator group_begin() const;
            /*!  \brief std::vector end iterator to recurse group names

                 \return End iterator
            */
            group_iterator group_end() const;

            /*!  \brief std::unordered_map iterator to recurse key-values

                 \return Start iterator
            */
            iterator begin();
            /*!  \brief std::unordered_map end iterator to recurse key-values

                 \return End iterator
            */
            iterator end();

            /*!  \brief std::unordered_map iterator to recurse key-values

                 \return Start iterator
            */
            const_iterator cbegin() const;
            /*!  \brief std::unordered_map end iterator to recurse key-values

                 \return End iterator
            */
            const_iterator cend() const;

            /*!  \brief focs::FilterIterator iterator to recurse key-values of the input group

                 \return Start iterator
            */
            group_filtered_iterator begin(const std::string& group) const;
            /*!  \brief focs::FilterIterator end iterator to recurse key-values of the input group

                 \return End iterator
            */
            group_filtered_iterator end(const std::string& group) const;


            /*!
                \brief If k matches the key of an element in the KvStore, the function returns a reference to its mapped value.

                If k does not match the key of any element in the container,
                the function inserts a new element with that key and returns a
                reference to its mapped value. Notice that this always
                increases the container size by one, even if no mapped value is
                assigned to the element (the element is constructed using its
                default constructor).

                \param k Key to look up

                \return A reference to the mapped value
            */
            std::string& operator[] (const std::string& k);
            /*! \copydoc operator[](const std::string&) */
            std::string& operator[] (std::string&& k);

            /*! \brief Returns a reference to the mapped value of the element identified with key k. */
            std::string& at(const std::string& k);
            /*! \copydoc at(const std::string&) */
            std::string& at(std::string&& k);

            /*! \brief Get iterator to element 

                Searches the container for an element with k as key and returns
                an iterator to it if found, otherwise it returns an iterator to
                unordered_map::end (the element past the end of the container).

                \param k Key to look up

                \return An iterator to the element, if the specified key value
                    is found, or unordered_map::end if the specified key is not
                    found in the container.
            */
            KvStore::iterator find(const std::string& k);
            /*! \copydoc find(const std::string&) */
            KvStore::iterator find(std::string&& k);
            /*! \copydoc find(const std::string&) */
            KvStore::const_iterator find(std::string& k) const;
            /*! \copydoc find(const std::string&) */
            KvStore::const_iterator find(std::string&& k) const;

            /*! \brief Count elements with a specific key

                Searches the container for an element with k as key and returns
                an iterator to it if found, otherwise it returns an iterator to
                unordered_map::end (the element past the end of the container).

                \param k Key to look up

                \return 1 if an element with a key equivalent to k is found, or zero otherwise.
            */
            size_t count(const std::string& k) const;
            /*! \copydoc count(const std::string&) */
            size_t count(std::string&& k) const;


            /*!
                \brief Output KvStore to stream operator

                \param os Target output stream
                \param kv KvStore to output

                \return Output stream, for chaining
            */
            friend std::ostream& operator<<(std::ostream& os, const KvStore& kv);

            /*!
                \brief Process a key-value store command by name

                \param cmd Command name to run
                \param args Arguments to pass to command function

                \return true on error, false on success
            */
            bool command(std::string&& cmd, std::string&& args);
            /*!
                \brief Change the current key group

                \param group New group name

                \return true on error, false on success
            */
            bool command_group(std::string&& group);
            /*!
                \brief Include another key-value store, relative to current file, into the global group

                \param file file relative to current key-value file path

                \return true on error, false on success
            */
            bool command_include(std::string&& file);
            /*!
                \brief Include another key-value store, relative to current file, into the current group

                \param file file relative to current key-value file path

                \return true on error, false on success
            */
            bool command_include_local(std::string&& file);

            /*!
                \brief Set the character(s) to use between group and key when inserting into map

                \param s New separator
            */
            void group_separator(const std::string& s){group_separator_ = s;}
            /*!
                \brief Get the current group separator

                \return the current group separator
            */
            const std::string& group_separator(){return group_separator_;}

            /*!
                \brief Set a function to be called for each key-value read

                \param callback New callback
            */
            void callback(read_callback callback){callback_ = callback;}
            /*!
                \brief Get the key-value callback

                \return the current key-value callback
            */
            read_callback callback(){return callback_;}

            /*!
                \brief Switch default group for retrieving values

                \param group Group to use by default

            */
            void switch_group(const std::string& group);

        private:
            std::vector<std::string> groups_{};
            std::unordered_map<std::string, std::string> kv_{};
            std::string group_separator_{"."};
            std::string current_group_{""};
            std::string current_file_{""};
            std::string read_group_{""};

            read_callback callback_{};

            /*!  \brief Add a group if it doesn't already exist

                 \param group new group
            */
            void add_group(const std::string& group);

    };
} // namespace focs

#endif // FOCS_KVSTORE

