// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.  Place custom code in custominit.h.
void mcExitClass_XTek_recon();
void mcInitClass_XTek_recon();
void mcExitClass_Parallel_Beam_recon();
void mcInitClass_Parallel_Beam_recon();
void mcExitClass_NeXus_normalise();
void mcInitClass_NeXus_normalise();
void mcExitClass_Cone_Beam_recon();
void mcInitClass_Cone_Beam_recon();


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void amirapackage_CCPiReconstruction_init()
{
    static bool isInitialized = false;
    if (isInitialized)
      return;

    isInitialized = true;

    mcInitClass_XTek_recon();
    mcInitClass_Parallel_Beam_recon();
    mcInitClass_NeXus_normalise();
    mcInitClass_Cone_Beam_recon();

}


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void amirapackage_CCPiReconstruction_finish()
{
    static bool isFinished = false;
    if (isFinished)
      return;

    isFinished = true;


    mcExitClass_Cone_Beam_recon();
    mcExitClass_NeXus_normalise();
    mcExitClass_Parallel_Beam_recon();
    mcExitClass_XTek_recon();
}

#if defined(_WIN32)
#  include <windows.h>
BOOL WINAPI DllMain(
    __in  HINSTANCE hinstDLL,
    __in  DWORD fdwReason,
    __in  LPVOID lpvReserved
    )
{
    switch (fdwReason)
    {
    case DLL_PROCESS_ATTACH:
        amirapackage_CCPiReconstruction_init();
        break;
    case DLL_PROCESS_DETACH:
        amirapackage_CCPiReconstruction_finish();
        break;
    default:
        ;
    }
    return true;
}


#endif

#if defined(__GNUC__)
void __attribute__((constructor)) soconstructor_CCPiReconstruction() {
    amirapackage_CCPiReconstruction_init();
}

void __attribute__((destructor)) sodestructor_CCPiReconstruction() {
    amirapackage_CCPiReconstruction_finish();
}
#endif
