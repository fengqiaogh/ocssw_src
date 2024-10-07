
#ifndef FOCS_MODULE
#define FOCS_MODULE

#include "DataProvider.hpp"
#include "FileReader.hpp"
#include "FileWriter.hpp"
#include "Log.hpp"

#include <filesystem>

#include <iostream>
#include <memory>
#include <string>
#include <vector>


// This is everything needed to create a new module

namespace focs {
    /*! The function in modules called to instantiate them.

        A function of this name must be found, exportable, with C linkage, that
        returns a newly allocated pointer to an object subclassing the
        focs::Module class.  Its only argument is a const reference to the
        ModuleConfiguration object.  Below is the definition used in the core
        module.

        focs::CoreModule* focs_module_allocator(const focs::ModuleConfiguration& configuration) {
            return new focs::CoreModule(configuration);
        }
    */
    static const std::string module_allocator_function = "focs_module_allocator";
    /*! The function in modules called to instantiate them.

        A function of this name must be found, exportable, with C linkage, that
        deletes the pointer returned from the allocator described above.  Below
        is the definition used in the core module.

        void focs_module_deleter(focs::CoreModule *module) {
            delete module;
        }
    */
    static const std::string module_deleter_function = "focs_module_deleter";

    class ModuleManager;
    /*! This configuration is passed from the master oc command to each module as it is is created. */
    struct ModuleConfiguration {
        /*! Sole constructor to initialize the object */
        ModuleConfiguration(ModuleManager& module_manager_, Log& log_) : module_manager{module_manager_}, log{log_} {}
        /*! An object containing references to all modules found (may not be full while constructing modules). */
        ModuleManager& module_manager;
        /*! The common logger that all modules should use and pass around */
        Log& log;
    };

    // TODO: Should this class get moved into its own header, like FileReader?
    // If helper classes/functions are needed, it will be.  (This will likely
    // happen, as it would behoove me to tie together KvStore and
    // getopt/boost::program_options/etc, or otherwise switch to JSON for
    // configs with extra processing for including other files, etc, like
    // KvStore has.  Must ask someone smarter.)

    /*! Superclass of commands to be added to oc */
    class Command {
    public:
        /*! Standard, empty virtual destructor */
        virtual ~Command() = default;
        /*! The name of the command, with optional aliases */
        virtual const std::vector<std::string>& respond_to() const = 0;
        /*! The function called when the command is run */
        virtual int call_command(int argc, const char* argv[]) const = 0;
        /*! Optional brief text summary, for documentation/help commands */
        virtual const std::string* brief_summary() const {return nullptr;}
        /*! Optional command group, for documentation/help commands */
        virtual const std::string* group() const {return nullptr;}
    };

    /*! Superclass of dynamically loaded modules */
    class Module {
    public:
        /*! Standard, empty virtual destructor */
        virtual ~Module() = default;

        /*! The name of the module, used only for display purposes */
        virtual const std::string& name() const = 0;
        /*! Optional, a list of commands the module provides, if any */
        virtual const std::vector<std::unique_ptr<Command>>& commands() const {return no_commands_;}
        /*! Optional, a list of data providers (AKA, algorithms) the module provides, if any */
        virtual const std::vector<std::unique_ptr<DataProvider>>& data_providers() const {return no_data_providers_;}
        /*! Optional, a list of file readers the module provides, if any */
        virtual const std::vector<std::unique_ptr<FileReader>>& file_readers() const {return no_file_readers_;}
        /*! Optional, a list of file writers the module provides, if any */
        virtual const std::vector<std::unique_ptr<FileWriter>>& file_writers() const {return no_file_writers_;}
    private:
        /*! Empty vector, for ease of default non-functionality */
        std::vector<std::unique_ptr<Command>> no_commands_{};
        /*! Empty vector, for ease of default non-functionality */
        std::vector<std::unique_ptr<DataProvider>> no_data_providers_{};
        /*! Empty vector, for ease of default non-functionality */
        std::vector<std::unique_ptr<FileReader>> no_file_readers_{};
        /*! Empty vector, for ease of default non-functionality */
        std::vector<std::unique_ptr<FileWriter>> no_file_writers_{};
    };
}

// Nothing below here is for user consumption.

namespace focs {
    class ModuleLoader {
    public:
        static bool is_module(const std::filesystem::path& path);

        ModuleLoader(const std::filesystem::path& path, const ModuleConfiguration& configuration);
        ~ModuleLoader();

        void load_module();
        std::shared_ptr<Module> get_module() const;
        void unload_module();
        bool is_loaded() const;
        const std::filesystem::path& path() const;

    private: // these need to accommodate every implementation or this needs to be split out into a base class again
        std::filesystem::path path_;
        const ModuleConfiguration& configuration_;
        void* handle_ { nullptr };
        std::shared_ptr<Module> module_;
    };
}

#endif // FOCS_MODULE

