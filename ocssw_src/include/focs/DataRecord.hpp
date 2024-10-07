
#ifndef FOCS_DATARECORD
#define FOCS_DATARECORD

#include "Variable.hpp"
#include "Product.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace focs {
struct TileParameters {
    public:
        TileParameters(size_t lines, size_t pixels, size_t start_line=0, size_t start_pixel=0) : size{lines, pixels}, origin{start_line, start_pixel}, end{origin.first + size.first, origin.second + size.second} {}
        TileParameters(std::pair<size_t, size_t> size_, std::pair<size_t, size_t> origin_={0,0}) : size{size_}, origin{origin_}, end{origin.first + size.first, origin.second + size.second} {}
        TileParameters(std::pair<size_t, size_t> size_, std::pair<size_t, size_t> origin_, std::pair<size_t, size_t> data_origin_, bool top_is_bounded_, bool bottom_is_bounded_, bool left_is_bounded_=true, bool right_is_bounded_=true) : 
            size{size_}, origin{origin_}, end{origin.first + size.first, origin.second + size.second}, data_origin{data_origin_}, top_is_bounded{top_is_bounded_}, bottom_is_bounded{bottom_is_bounded_}, left_is_bounded{left_is_bounded_}, right_is_bounded{right_is_bounded_} {}

        void shrink_to(std::pair<size_t, size_t> new_size){
            size.first = std::min(size.first, new_size.first);
            size.second = std::min(size.second, new_size.second);
            end = {origin.first + size.first, origin.second + size.second};
        }
        void offset_up_to(std::pair<size_t, size_t> new_position){
            origin.first = std::max(origin.first, new_position.first);
            origin.second = std::max(origin.second, new_position.second);
            end = {origin.first + size.first, origin.second + size.second};
        }
        void update(const TileParameters& params){
            shrink_to(params.size);
            offset_up_to(params.origin);
            data_origin = params.data_origin;
            top_is_bounded |= params.top_is_bounded;
            bottom_is_bounded |= params.bottom_is_bounded;
            left_is_bounded |= params.left_is_bounded;
            right_is_bounded |= params.right_is_bounded;
        }
        bool is_empty(){
            return size.first == 0 || size.second == 0;
        }

        friend std::ostream& operator<<(std::ostream& out, const TileParameters& params){
            out << "TileParameters[" << params.size.first << "x" << params.size.second << ";" << params.origin.first << "," << params.origin.second << ";";
            if (params.top_is_bounded){
                out << "^";
            }
            if (params.bottom_is_bounded){
                out << "v";
            }
            if (params.left_is_bounded){
                out << "<";
            }
            if (params.right_is_bounded){
                out << ">";
            }
            out << ";" << params.data_origin.first << "," << params.data_origin.second;
            return out << "]";
        }

        std::pair<size_t, size_t> size{0,0};
        std::pair<size_t, size_t> origin{0,0};
        std::pair<size_t, size_t> end{0,0}; // quality-of-life, origin + size
        std::pair<size_t, size_t> data_origin{0,0};
        bool top_is_bounded{false};
        bool bottom_is_bounded{false};
        bool left_is_bounded{false};
        bool right_is_bounded{false};
};

// unused
class DataRecordProcessor {
    public:
    private:
        std::unique_ptr<BaseVariable> input_{};
        std::unique_ptr<BaseVariable> output_{};
};

class DataRecord {
    public:
        DataRecord(){}
        // BaseVariable* get_product(const std::string& product) { return variables_.at(product).get(); }

        std::unordered_map<std::string, const Product>::const_iterator products_begin() const { return products_.cbegin(); }
        std::unordered_map<std::string, const Product>::const_iterator products_end() const { return products_.cend(); }

        std::pair<size_t, size_t> size() const noexcept { return size_; }
        void size(std::pair<size_t, size_t> size) noexcept { size_ = size; }
        std::pair<size_t, size_t> max_size() const noexcept { return max_size_; }
        void max_size(std::pair<size_t, size_t> max_size) noexcept { max_size_ = max_size; }
        void add(std::unique_ptr<BaseVariable> v){
            variables_.push_back(std::move(v));
        }
        auto& attributes() { return attributes_; }
        auto& dimensions() { return dimensions_; }
        auto& variables() { return variables_; }
        void min_size(size_t min_size){ min_size_ = std::max(min_size, min_size_); }
        size_t min_size() const { return min_size_; }

        // Reserves to record size, not data dimensions
        template<typename T>
        void reserve_data(std::vector<std::vector<T>>& data, const bool resize=true) const {
            data.reserve(size_.first);
            for (size_t i=0;i<size_.first;i++){
                data.emplace_back(0);
            }
            if (resize){
                for (auto& v : data){
                    v.resize(size_.second);
                }
            } else {
                for (auto& v : data){
                    v.reserve(size_.second);
                }
            }
        }

        // TODO: should the create_variable functions be extended for interpolation and/or extrapolation?
        // They'll accept two sets of dimensions and a conversion method.
        // They will create a hidden variable of the source dimensions and return a real one based on the second set of dimensions.
        // During the done_reading_input method, the conversion will be performed.
        // 1. In what form is the conversion method?
        //      enum and 100% controlled here?
        //      object?  More work for users, less defaulting.
        // 2. How are 2D interps done?  The kernel size is going to effect most
        //      maths.  It'll be up to the input file to set it, but a much larger
        //      min_size_ could be chosen.  8 or 16 rows at a time should not be too
        //      memory heavy, right?
        template<typename T, size_t N>
        focs::Variable<T, N>* create_variable(focs::Product& product, const std::vector<std::pair<std::string, size_t>>& dimensions={}, const bool resize=true) {
            return create_variable<focs::Variable<T, N>>(product, dimensions, resize);
        }

        template<typename T=focs::Variable<float, 2>>
        T* create_variable(focs::Product& product, const std::vector<std::pair<std::string, size_t>>& dimensions={}, const bool resize=true) {
            auto new_var = std::make_unique<T>(product);
            reserve_data(new_var->data(), resize);

            for (const auto& dim : dimensions){
                new_var->dimensions().emplace_back(dim.first, dim.second);
            }

            T* ret = new_var.get();
            product.variable(ret);
            variables_.push_back(std::move(new_var));
            return ret;
        }

        void done_reading_input(){
        }

        void rotate_data(){
            if (min_size() > 1){
                const size_t rotation_size = (min_size() >> 1) * 2;
                // std::cout << "Rotating " << rotation_size << " lines\n";
                for (auto& v : variables_){
                    // std::cout << "Rotating variable " << v->name() << "\n";
                    v->rotate(rotation_size);
                }
            }
        }

    private:
        std::unordered_map<std::string, size_t> dimensions_{};
        // std::vector<DataRecordProcessor> variables_to_process_{};
        std::vector<std::unique_ptr<BaseVariable>> variables_{};
        std::unordered_map<std::string, const Product> products_{};
        std::vector<std::unordered_map<std::string, const Attribute>> attributes_{}; // file attributes?

        std::pair<size_t, size_t> size_{0,0};
        std::pair<size_t, size_t> max_size_{0,0};
        size_t min_size_{1};
};

} // namespace focs

#endif // FOCS_DATARECORD
