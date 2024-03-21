#ifndef DATA_STORAGE_HPP
#define DATA_STORAGE_HPP

#include <hdf5.h>
#include <string>
#include <type_traits>
#include <vector>

// Some definitions for convenience
template <typename T> using Invoke = typename T::type;
template <typename Condition> using negation = std::integral_constant<bool, !bool(Condition::value)>;
template <typename Condition> using EnableIf = Invoke<std::enable_if<Condition::value>>;
template <typename Condition> using DisableIf = EnableIf<negation<Condition>>;

// Do not deduce from templates for pointers
#define template_noPointer template <typename T, typename = DisableIf<std::is_pointer<T>>>

template <typename T> static hid_t get_hdf5_data_type();

template <> hid_t get_hdf5_data_type<int>() { return H5Tcopy(H5T_NATIVE_INT); }
template <> hid_t get_hdf5_data_type<float>() { return H5Tcopy(H5T_NATIVE_FLOAT); }
template <> hid_t get_hdf5_data_type<double>() { return H5Tcopy(H5T_NATIVE_DOUBLE); }

class data_storage {
public:
	data_storage(std::string file_name);
	void write_dataset(const std::vector<double> &raw_data, const std::vector<size_t> &extents, std::string dataset_name);
	template_noPointer void AddGlobalAttr(const std::string &AttrName, T AttrData);

private:
	void open_file();
	void close_file();
	hid_t hdf5file, hdf5group;
	std::string file_name;
};

#endif