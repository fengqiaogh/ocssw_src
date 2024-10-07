
#ifndef FOCS_FILEREADER
#define FOCS_FILEREADER

#include <array>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "optional"
#include "DataProvider.hpp"
#include "DataRecord.hpp"
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/ostreamwrapper.h>
#include "rapidjson/error/en.h"
#include <netcdf>
#include "productInfo.h"
#include "sensorInfo.h"
#include <stack>
#include <l1.h>
namespace focs {
namespace rj = rapidjson;
namespace fs = std::filesystem;
extern const std::array<int, 4> NETCDF1_MAGIC_NUMBER;
extern const std::array<int, 4> NETCDF2_MAGIC_NUMBER;
extern const std::array<int, 4> HDF4_MAGIC_NUMBER;
extern const std::array<int, 8> HDF5_MAGIC_NUMBER;

struct ProdAttrGeneric {
    std::unordered_map<std::string, std::unordered_map<std::string, ANY>> attributes;
    std::set<std::string> product_names;
};

template <typename, typename = void>
constexpr bool is_prod_gen_attr{};

template <typename T>
constexpr bool
    is_prod_gen_attr<T, std::void_t<decltype(std::declval<T>().getAttrGen())>> =
        true;

class ConfigReader {
   protected:
    std::unique_ptr<rj::Document> document;
    std::optional<fs::path> path;
    [[maybe_unused]] std::unique_ptr<productInfo_t> product_xml;  // from ocssw oel_util, not implemented yet
    ProdAttrGeneric attrs;

   public:
    ConfigReader() = default;
    ConfigReader(const std::filesystem::path& path);
    ConfigReader(const fs::path& path,const std::string& type); // generic json or xml file
    bool exist() const;
    const rj::Document* getDoc() const;
    const ProdAttrGeneric & getAttrGen() const;
};

// Everything that writes to hints.attributes_ is performing a copy (of what they're adding).  Switch to smart pointers?
class FileReaderHints {
    public:
        FileReaderHints() = default;
        explicit FileReaderHints(const std::string& path);

        bool magic_number(const std::string& bytes) const noexcept;
        template<std::size_t N> bool magic_number(const std::array<int, N>& bytes) const noexcept;

        const std::array<int, 8>& first_bytes() const noexcept {return first_bytes_;}
        void first_bytes(const std::array<int, 8>& bytes){first_bytes_ = bytes;}

        uintmax_t size() const noexcept {return size_;}
        void size(uintmax_t size){size_ = size;}

        const std::string& filename() const noexcept {return filename_;}
        void filename(const std::string& filename){filename_ = filename;}

        const std::string& format() const noexcept {return format_;}
        void format(const std::string& format){format_ = format;}

        auto& attributes(){return attributes_;}

        bool is_hdf4() const noexcept;
        bool is_hdf5() const noexcept;
        bool is_netcdf() const noexcept;
        std::set<Attribute> parse_hdf4_metadata(const std::string &);
        // hdf4 legacy keywords
        const static std::string  core_metadata;
        const static std::string  struct_metadata;
        const static std::string  archive_metadata;
    private:
        std::string filename_{};
        std::array<int, 8> first_bytes_{{0, 0, 0, 0, 0, 0, 0, 0}};
        uintmax_t size_{0};
        std::string format_{"unknown"};
        std::set<focs::Attribute> attributes_{};

         
};
class FileReader : public DataProvider {
    public:
        FileReader(const std::string& name, const std::string& description) : DataProvider(name, description) {}
        FileReader(const std::string& name) : FileReader(name, {}) {}

        enum validity {
            invalid, possibly_valid, valid
        };
        virtual ~FileReader() override {}

        // virtual const std::vector<Product>& provides() const override;
        virtual std::vector<std::shared_ptr<Product>>& needs() override; // returns empty

        // TODO: undefine these empty ones, it's not compiling without them for some reason
        virtual validity is_valid_file(const std::string& file, FileReaderHints& hints){(void)file; (void)hints; return invalid;}
        validity is_valid_file(const std::string& file){
            FileReaderHints hints{file};
            return is_valid_file(file, hints);
        }

        virtual std::unique_ptr<FileReader> initialize_reader(DataProviderConfiguration& configuration, const std::filesystem::path& path) = 0;
        virtual std::set<focs::Attribute> get_file_attributes() = 0;
        virtual TileParameters read_next_tile(DataProviderConfiguration& configuration, DataRecord& record){(void)configuration; (void)record; return {0,0}; } // true on finished? For now.
};


class AttributeReader 
{
protected:
    std::unique_ptr<ConfigReader> _config;
    std::unordered_map<std::string,ConfigReader> configs;
    std::unordered_map<std::string,std::set<focs::Attribute>> _attributes;
    productInfo_t *info = nullptr;

   public:
    AttributeReader(const std::filesystem::path& path, const std::string& mode = "default",
                    const std::unordered_set<std::string>& products_requested = {});
    AttributeReader(netCDF::NcFile& nc, const std::string& key);
    AttributeReader(FileReaderHints& hints, const std::string& key);
    AttributeReader();
    ~AttributeReader() = default;
    void addAttributes(const std::filesystem::path& path, const std::string& mode = "default",
                       const std::unordered_set<std::string>& products_requested = {});
    void addAttributes(FileReaderHints& hints, const std::string& key);
    void addAttributes(netCDF::NcFile& nc, const std::string& key);
    template <typename T, std::enable_if_t<is_prod_gen_attr<T>, bool> = true,
              typename pointer_type = decltype((std::declval<T>().getAttrGen())),
              std::enable_if_t<std::is_same_v<pointer_type, const ProdAttrGeneric&>, bool> = true>
    void addAttributes(T&& storage_object, const std::string& mode,
                       const std::unordered_set<std::string>& products_requested = {}) {
        const ProdAttrGeneric& prods = storage_object.getAttrGen();
        if (mode == "add") {
            for (const auto& prod_name : prods.product_names) {
                if (!products_requested.empty() && products_requested.count(prod_name) == 0)
                    continue;
                if (_attributes.count(prod_name) > 0) {
                    std::unordered_set<std::string> existing_attrs;
                    for (const auto& attr : _attributes.at(prod_name)) {
                        existing_attrs.insert(attr.name());
                    }
                    // insert ANY attributes
                    for (const auto& attr : prods.attributes.at(prod_name)) {
                        std::string name = attr.first;
                        if (existing_attrs.count(name) == 0) {
                            _attributes.at(prod_name).insert(focs::Attribute{name, attr.second, 1});
                        }
                    }
                } else {
                    for (const auto& attr : prods.attributes.at(prod_name)) {
                        _attributes[prod_name].insert(focs::Attribute{attr.first, attr.second, 1});
                    }
                }
            }
        } else if (mode == "override") {
            // needs to be implemented
        }
    }
    template <typename T, std::enable_if_t<is_iterable<T>, bool> = true,
              typename pointer_type = decltype(*(std::declval<T>().begin())),
              std::enable_if_t<std::is_same_v<pointer_type, const Attribute&>, bool> = true>  //
    void addAttributes(const T& custom_attrs, const std::string& type, const std::string& mode = "add") {
        if (mode == "add") {
            std::unordered_set<std::string> existing_attrs;
            for (const auto& attr : _attributes.at(type)) {
                existing_attrs.insert(attr.name());
            }
            for (const auto& attr : custom_attrs) {
                std::string name = attr.name();
                if (existing_attrs.count(name) == 0) {
                    _attributes.at(type).insert(attr);
                }
            }

        } else if (mode == "override") {
            std::unordered_set<std::string> override_names;
            for (const auto& attr : custom_attrs) {
                override_names.insert(attr.name());
            }
            std::stack<std::set<focs::Attribute>::iterator> its_to_remove;
            for (auto it = _attributes.at(type).begin(); it != _attributes.at(type).end(); it++) {
                if (override_names.count(it->name()) > 0) {
                    its_to_remove.push(it);
                }
            }
            while (!its_to_remove.empty()) {
                auto it = its_to_remove.top();
                _attributes.at(type).erase(it);
                its_to_remove.pop();
            }
            _attributes.at(type).insert(custom_attrs.begin(), custom_attrs.end());
        }
    }
    const  std::unordered_map<std::string,std::set<focs::Attribute>> & getAttrs() const;
};

} // namespace focs

#endif // FOCS_FILEREADER

