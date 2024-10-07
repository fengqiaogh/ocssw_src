
#ifndef FOCS_MODULEMANAGER
#define FOCS_MODULEMANAGER

#include "Module.hpp"

#include <boost/filesystem.hpp>

#include <memory>
#include <vector>

namespace focs {
    class ModuleManager {
    public:
        ModuleManager(Log& log);
        ModuleManager(const std::filesystem::path& module_directory, Log& log);
        ~ModuleManager();

        void load_modules();
        const std::vector<std::shared_ptr<Module>>& get_modules() const;
        const std::vector<ModuleLoader>& get_loaders() const;

        static const std::filesystem::path default_module_directory;
    private:
        std::filesystem::path module_directory_{default_module_directory};
        std::vector<std::shared_ptr<Module>> modules_{};
        std::vector<ModuleLoader> module_loaders_{};
        // bool is_loaded_{false};
        Log& log_;
        ModuleConfiguration configuration_;
    };
}

#endif // FOCS_MODULEMANAGER
