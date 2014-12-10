/// \page cmd_3d_hybrid_reconstruction 3d_hybrid_reconstruction
/// \brief 3D Hybrid PET system reconstruction tool
///
/// NYI
/// ===

#include "detector_set.h"
#include "options.h"

using namespace PET3D::Hybrid;

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    add_reconstruction_options(cl);
    cl.try_parse(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    return 0;
  } catch (cmdline::exception& ex) {
    if (ex.help()) {
      std::cerr << ex.usage();
    }
    for (auto& msg : ex.errors()) {
      auto name = ex.name();
      if (name) {
        std::cerr << "error at " << name << ": " << msg << std::endl;
      } else {
        std::cerr << "error: " << msg << std::endl;
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}
