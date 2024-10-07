#ifndef FOCS_DATAPROVIDER
#define FOCS_DATAPROVIDER

#include "DataRecord.hpp"
#include "KvStore.hpp"
#include "Product.hpp"

#include <algorithm>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <iostream> 

namespace focs {
    class DataProvider; // circular dependency
    struct DataProviderPathPart {
        DataProvider* provider{nullptr};
        std::unordered_set<Product*> provide{};
        std::unordered_set<std::shared_ptr<Product>> need{};
        bool is_input{false}; // obsolete?
        bool is_output{false}; // obsolete?

        bool provides(const std::shared_ptr<Product> needle) const {
            return std::find_if(provide.cbegin(), provide.cend(), [&needle](const auto* haystack){return haystack->matches(*needle);}) != provide.cend();
        }
        void update_needs(std::vector<std::shared_ptr<Product>> new_needs) {
            for (auto& n : new_needs){
                add_or_find_need(n);
            }
        }
        std::shared_ptr<Product> add_or_find_need(std::shared_ptr<Product> new_need) {
            auto it = std::find_if(need.cbegin(), need.cend(), [&new_need](const auto& this_need){return *this_need == *new_need;});
            if (it == need.cend()){
                need.insert(new_need);
                return new_need;
            }
            return *it;
        }
        DataProviderPathPart(DataProvider* provider_) : provider{provider_} {}
        DataProviderPathPart(DataProvider* provider_, std::unordered_set<Product*> provide_) : provider{provider_}, provide{provide_} {
            std::cout << "Creating path part, " << provide.size() << "\n";
        }
        DataProviderPathPart(DataProvider* provider_, bool is_input_) : provider{provider_}, is_input{is_input_} {} // obsolete
        DataProviderPathPart(DataProvider* provider_, bool is_input_, bool is_output_) : provider{provider_}, is_input{is_input_}, is_output{is_output_} {} // obsolete
    };

    // given to initialize data providers after a processing path has been chosen
    class DataProviderConfiguration {
        public:
            DataProviderConfiguration(DataRecord& data_record_,
            std::vector<std::shared_ptr<DataProviderPathPart>>& processing_path_, std::unordered_map<std::string,std::string> & configuration_) : configuration{configuration_}, data_record{data_record_},processing_path{processing_path_} {}

            KvStore configuration{};
            DataRecord& data_record;
            std::vector<std::shared_ptr<DataProviderPathPart>>& processing_path;
    };
    class DataProvider {
        public:
            DataProvider(const std::string& name, const std::string& description) : name_{name}, description_{description} {}
            DataProvider(const std::string& name) : DataProvider(name, {}) {}

            virtual ~DataProvider();

            const std::string& name() const {return name_;}
            const std::string& description() const {return description_;}

            // virtual std::unique_ptr<DataProvider> create(DataProviderConfiguration& configuration) = 0;

            virtual void pre_initialize(DataProviderConfiguration& configuration){(void)configuration;}
            virtual void initialize(DataProviderConfiguration& configuration, std::unordered_set<Product*>& provide){(void)configuration;(void)provide;}
            virtual void post_initialize(DataProviderConfiguration& configuration){(void)configuration;}
            virtual void finalize(){}

            virtual void process_tile(DataProviderConfiguration& configuration, DataRecord &data_record, TileParameters& tile_params) = 0;

            virtual std::vector<Product>& provides();
            virtual void set_needs(std::unordered_set<std::shared_ptr<Product>>& provide){(void)provide;} // only providers that dynamically create needs/provides during path searching need to utilize this
            virtual std::vector<std::shared_ptr<Product>>& needs();
            virtual std::vector<std::shared_ptr<Product>> needs(std::unordered_set<Product*>& provide){ // TODO:  should this get renamed?
                (void)provide;
                return needs();
            }

            virtual size_t kernel_size() const { return 1; }


            //TODO: fulfills needs to be overloadable for providers that only need to check certain attributes.
            // These won't have to be a perfect match because they may be able to assume that the rest of the attributes will be filled in by whatever they require.

            virtual bool fulfills(const Product& other) {
                return contains(provides(), other);
            }
            bool fulfills(const std::vector<std::shared_ptr<Product>>& other) {
                return contains_all(provides(), other);
            }
            bool fulfills(DataProvider& other) {
                return contains_all(provides(), other.needs());
            }
            bool is_fulfilled_by(const std::vector<Product>& other) {
                return contains_all(other, needs());
            }
            bool is_fulfilled_by(DataProvider& other) {
                return contains_all(other.provides(), needs());
            }
            bool partially_fulfills(const Product& other) { // convenience
                return fulfills(other);
            }
            bool partially_fulfills(const std::vector<std::shared_ptr<Product>>& other) {
                return contains_any(provides(), other);
            }
            bool partially_fulfills(DataProvider& other) {
                return contains_any(provides(), other.needs());
            }
            bool partially_fulfills(const std::vector<DataProvider*> other) {
                return std::any_of(other.cbegin(), other.cend(), [&](const auto& o){ return contains_any(this->provides(), o->needs()); });
            }
            bool is_partially_fulfilled_by(const std::vector<Product>& other) {
                return contains_any(other, needs());
            }
            bool is_partially_fulfilled_by(DataProvider& other) {
                return contains_any(other.provides(), needs());
            }
            // This one is weird, see is_contained_in below
            bool is_partially_fulfilled_by(const Product& other) {
                return is_contained_in(other, needs());
            }

            static std::vector<std::shared_ptr<DataProviderPathPart>> processing_paths(
                    KvStore& configuration, 
                    std::vector<DataProvider*>& input, 
                    std::vector<DataProvider*>& data_providers, 
                    std::vector<DataProvider*>& output);

        private:
            std::string name_{};
            std::string description_{};

            static bool contains_all(const std::vector<Product>& haystack, const std::vector<std::shared_ptr<Product>>& needles) {
                for (auto& needle : needles){
                    if (!contains(haystack, *needle)){
                        return false;
                    }
                }
                return true;
            }
            static bool contains_any(const std::vector<Product>& haystack, const std::vector<std::shared_ptr<Product>>& needles) {
                for (auto& needle : needles){
                    if (contains(haystack, *needle)){
                        return true;
                    }
                }
                return false;
            }
            static bool contains(const std::vector<Product>& haystack, const Product& needle) {
                for (auto& hay : haystack){
                    if (hay.matches(needle)){
                        return true;
                    }
                }
                return false;
            }

            // This naming is a bit weird due to Product::match being
            // asymmetrical.  The only time this is used, we're checking if any
            // needed product is fulfilled by an input product that is
            // provided.
            static bool is_contained_in(const Product& hay, const std::vector<std::shared_ptr<Product>>& needles){
                for (auto& needle : needles){
                    if (hay.matches(*needle)){
                        return true;
                    }
                }
                return false;
            }
    };

} // namespace focs

std::ostream& operator<<(std::ostream& out, focs::DataProvider& in);

#endif // FOCS_DATAPROVIDER

