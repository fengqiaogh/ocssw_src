#ifndef FOCS_FILEWRITER
#define FOCS_FILEWRITER

#include <array>
#include <cstdio>
#include <string>
#include <unordered_map>

#include <boost/filesystem.hpp>

#include "DataProvider.hpp"
#include "DataRecord.hpp"

namespace focs {

class FileWriter : public DataProvider {
    public:
        FileWriter(const std::string& name, const std::string& description) : DataProvider(name, description) {}
        FileWriter(const std::string& name) : FileWriter(name, {}) {}

        virtual ~FileWriter() override {}

        virtual std::vector<Product>& provides() override; // returns empty
        // virtual std::vector<Product>& needs() override;


        virtual std::unique_ptr<FileWriter> initialize_writer(DataProviderConfiguration& configuration, const std::string& group) = 0;
        virtual void set_global_attributes(const std::set<focs::Attribute> &) {};
        virtual void set_product_attributes(const std::string&, const std::set<focs::Attribute> &) {};
        enum validity {
            invalid, possibly_valid, valid
        };
        virtual validity can_process_output_group(focs::DataProviderConfiguration& configuration, const std::string& group){
            (void)configuration; (void)group;
            return invalid;
        }

};


} // namespace focs

#endif // FOCS_FILEWRITER

