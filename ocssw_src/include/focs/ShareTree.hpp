#ifndef FOCS_SHARETREE
#define FOCS_SHARETREE

#include <filesystem>

#include <rapidjson/document.h>

#include <iterator>
#include <string>
#include <map>

namespace focs {
    struct SensorDirectory {
        std::string instrument;
        std::string platform;
        std::filesystem::path root;

        SensorDirectory(const rapidjson::GenericValue<rapidjson::UTF8<>>& in);
    };
    std::ostream& operator<<(std::ostream& out, const SensorDirectory& in);

    class ShareTree {
        public:
            ShareTree();
            ShareTree(std::filesystem::path root);

            inline static ShareTree& get_default(){
                static ShareTree _instance{};
                return _instance;
            }

            std::map<const std::string, const SensorDirectory>::const_iterator begin() const;
            std::map<const std::string, const SensorDirectory>::const_iterator end() const;

            const SensorDirectory& operator[] (const std::string& k) const;
            const SensorDirectory& operator[] (std::string&& k) const;

            const SensorDirectory& at(const std::string& k) const;

            std::string get_path_to_file(const std::string& sensor_id, const std::filesystem::path& rel_path);
            std::string get_path_to_file(const std::filesystem::path& rel_path);
        private:
            const std::filesystem::path root_;
            const std::map<const std::string, const SensorDirectory> directories_;
    };
} // namespace focs

#endif // FOCS_SHARETREE

