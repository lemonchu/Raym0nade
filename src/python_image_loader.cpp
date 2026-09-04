#include "python_image_loader.hpp"

#include <Python.h>

#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>

namespace raym0nade::detail {
namespace {

struct PythonObjectDeleter {
    void operator()(PyObject* object) const noexcept {
        Py_XDECREF(object);
    }
};

using PythonObject = std::unique_ptr<PyObject, PythonObjectDeleter>;

[[noreturn]] void throwPythonError(const std::string& message) {
    if (PyErr_Occurred() != nullptr) {
        PyErr_Print();
    }
    throw std::runtime_error(message);
}

class GilGuard {
public:
    GilGuard() : state_(PyGILState_Ensure()) {}

    GilGuard(const GilGuard&) = delete;
    GilGuard& operator=(const GilGuard&) = delete;

    ~GilGuard() {
        PyGILState_Release(state_);
    }

private:
    PyGILState_STATE state_;
};

class BufferView {
public:
    explicit BufferView(PyObject* object, int flags) {
        if (PyObject_GetBuffer(object, &view_, flags) != 0) {
            throwPythonError("The Python decoder returned a non-contiguous image buffer.");
        }
        acquired_ = true;
    }

    BufferView(const BufferView&) = delete;
    BufferView& operator=(const BufferView&) = delete;

    ~BufferView() {
        if (acquired_) {
            PyBuffer_Release(&view_);
        }
    }

    [[nodiscard]] const Py_buffer& get() const noexcept {
        return view_;
    }

private:
    Py_buffer view_{};
    bool acquired_{false};
};

void registerScriptDirectory() {
#ifdef RAYM0NADE_SCRIPT_DIR
    const std::filesystem::path scriptDirectory{RAYM0NADE_SCRIPT_DIR};
#else
    const std::filesystem::path scriptDirectory = std::filesystem::current_path() / "scripts";
#endif

    PyObject* systemPath = PySys_GetObject("path");
    if (systemPath == nullptr) {
        throwPythonError("Failed to access Python's module search path.");
    }

    const std::string encodedPath = scriptDirectory.u8string();
    if (encodedPath.size() > static_cast<std::size_t>(PY_SSIZE_T_MAX)) {
        throw std::runtime_error("The image-decoder script path is too long for Python.");
    }
    PythonObject pythonPath{PyUnicode_DecodeUTF8(
        encodedPath.data(), static_cast<Py_ssize_t>(encodedPath.size()), nullptr)};
    if (!pythonPath || PyList_Insert(systemPath, 0, pythonPath.get()) != 0) {
        throwPythonError("Failed to register the image-decoder script directory.");
    }
}

void initializePython() {
    static std::once_flag initializationFlag;
    std::call_once(initializationFlag, [] {
        const bool createdInterpreter = Py_IsInitialized() == 0;
        if (createdInterpreter) {
            Py_Initialize();
        }
        if (Py_IsInitialized() == 0) {
            throw std::runtime_error("Failed to initialize the embedded Python interpreter.");
        }

        if (createdInterpreter) {
            try {
                registerScriptDirectory();
            } catch (...) {
                (void)PyEval_SaveThread();
                throw;
            }
            (void)PyEval_SaveThread();
        } else {
            GilGuard guard;
            registerScriptDirectory();
        }
    });
}

PythonObject invokeDecoder(const std::filesystem::path& filename, const std::string& moduleName) {
    PythonObject pythonModuleName{PyUnicode_FromString(moduleName.c_str())};
    if (!pythonModuleName) {
        throwPythonError("Failed to encode the Python decoder module name.");
    }

    PythonObject module{PyImport_Import(pythonModuleName.get())};
    if (!module) {
        throwPythonError("Failed to import Python image decoder: " + moduleName);
    }

    PythonObject function{PyObject_GetAttrString(module.get(), moduleName.c_str())};
    if (!function || PyCallable_Check(function.get()) == 0) {
        throwPythonError("Python image decoder is not callable: " + moduleName);
    }

    const std::string encodedFilename = filename.u8string();
    if (encodedFilename.size() > static_cast<std::size_t>(PY_SSIZE_T_MAX)) {
        throw std::runtime_error("The image filename is too long for Python.");
    }
    PythonObject pythonFilename{PyUnicode_DecodeUTF8(
        encodedFilename.data(), static_cast<Py_ssize_t>(encodedFilename.size()), nullptr)};
    if (!pythonFilename) {
        throwPythonError("Failed to encode the image filename for Python.");
    }

    PythonObject arguments{PyTuple_Pack(1, pythonFilename.get())};
    if (!arguments) {
        throwPythonError("Failed to create Python decoder arguments.");
    }

    PythonObject result{PyObject_CallObject(function.get(), arguments.get())};
    if (!result) {
        throwPythonError("Python image decoder failed for: " + filename.string());
    }
    if (PyTuple_Check(result.get()) == 0 || PyTuple_Size(result.get()) != 4) {
        throw std::runtime_error("Python image decoder returned an invalid result tuple.");
    }
    return result;
}

int tupleInteger(PyObject* tuple, Py_ssize_t index, const char* fieldName) {
    PyObject* value = PyTuple_GetItem(tuple, index);
    if (value == nullptr || value == Py_None || PyLong_Check(value) == 0) {
        throw std::runtime_error(std::string("Python image decoder returned an invalid ") + fieldName + '.');
    }

    const long converted = PyLong_AsLong(value);
    if (converted <= 0 || converted > std::numeric_limits<int>::max() || PyErr_Occurred() != nullptr) {
        throwPythonError(std::string("Python image decoder returned an out-of-range ") + fieldName + '.');
    }
    return static_cast<int>(converted);
}

std::size_t checkedElementCount(int width, int height, int channels) {
    if (channels < 1 || channels > 4) {
        throw std::runtime_error("Image decoder returned an unsupported channel count.");
    }

    const auto imageWidth = static_cast<std::size_t>(width);
    const auto imageHeight = static_cast<std::size_t>(height);
    if (imageWidth > std::numeric_limits<std::size_t>::max() / imageHeight) {
        throw std::runtime_error("Decoded image dimensions overflow addressable memory.");
    }
    const auto pixelCount = imageWidth * imageHeight;
    if (pixelCount > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(channels)) {
        throw std::runtime_error("Decoded image dimensions overflow addressable memory.");
    }
    return pixelCount * static_cast<std::size_t>(channels);
}

}  // namespace

ByteImage loadByteImage(const std::filesystem::path& filename, const std::string& decoderModule) {
    initializePython();
    GilGuard guard;
    PythonObject result = invokeDecoder(filename, decoderModule);

    ByteImage image;
    image.width = tupleInteger(result.get(), 0, "width");
    image.height = tupleInteger(result.get(), 1, "height");
    image.channels = tupleInteger(result.get(), 2, "channel count");
    const std::size_t expectedSize = checkedElementCount(image.width, image.height, image.channels);

    PyObject* data = PyTuple_GetItem(result.get(), 3);
    if (data == nullptr || data == Py_None) {
        throw std::runtime_error("Python image decoder returned no pixel data.");
    }

    BufferView buffer{data, PyBUF_C_CONTIGUOUS};
    if (buffer.get().buf == nullptr || buffer.get().len < 0 ||
        static_cast<std::size_t>(buffer.get().len) != expectedSize) {
        throw std::runtime_error("Decoded byte image buffer size does not match its dimensions.");
    }

    image.pixels.resize(expectedSize);
    std::memcpy(image.pixels.data(), buffer.get().buf, expectedSize);
    return image;
}

FloatImage loadFloatImage(const std::filesystem::path& filename, const std::string& decoderModule) {
    initializePython();
    GilGuard guard;
    PythonObject result = invokeDecoder(filename, decoderModule);

    FloatImage image;
    image.width = tupleInteger(result.get(), 0, "width");
    image.height = tupleInteger(result.get(), 1, "height");
    image.channels = tupleInteger(result.get(), 2, "channel count");
    const std::size_t expectedSize = checkedElementCount(image.width, image.height, image.channels);
    if (expectedSize > std::numeric_limits<std::size_t>::max() / sizeof(float)) {
        throw std::runtime_error("Decoded HDR dimensions overflow addressable memory.");
    }

    PyObject* data = PyTuple_GetItem(result.get(), 3);
    if (data == nullptr || data == Py_None) {
        throw std::runtime_error("Python image decoder returned no pixel data.");
    }

    BufferView buffer{data, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT};
    const std::size_t expectedBytes = expectedSize * sizeof(float);
    const bool isFloat32 = buffer.get().itemsize == static_cast<Py_ssize_t>(sizeof(float)) &&
                           buffer.get().format != nullptr && std::strcmp(buffer.get().format, "f") == 0;
    if (!isFloat32 || buffer.get().buf == nullptr || buffer.get().len < 0 ||
        static_cast<std::size_t>(buffer.get().len) != expectedBytes) {
        throw std::runtime_error("Decoded HDR buffer is not a contiguous float32 image.");
    }

    image.pixels.resize(expectedSize);
    std::memcpy(image.pixels.data(), buffer.get().buf, expectedBytes);
    return image;
}

}  // namespace raym0nade::detail
