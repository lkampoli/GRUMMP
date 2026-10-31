// tri_libs.cpp (tri.exe)

//---------------------------------------------------------
// link with GRUMMP 2D libraries (v0.3.0)
//---------------------------------------------------------

#ifdef _DEBUG
# pragma message("Linking with GRUMMP v0.3.0 libraries (2D) -- static DEBUG build ")
# pragma comment(lib,"GR_baseD.lib")
# pragma comment(lib,"GR_2DD.lib")
# pragma comment(lib,"GR_3DD.lib")
# pragma comment(lib,"OptMSD.lib")
# pragma comment(lib,"SUMAAlogD.lib")
#else
# pragma message("Linking with GRUMMP v0.3.0 libraries (2D) -- static RELEASE build ")
# pragma comment(lib,"GR_base.lib")
# pragma comment(lib,"GR_2D.lib")
# pragma comment(lib,"GR_3D.lib")
# pragma comment(lib,"OptMS.lib")
# pragma comment(lib,"SUMAAlog.lib")
#endif
