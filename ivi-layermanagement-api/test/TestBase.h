
#include "wayland-client.h"
#include <vector>

extern "C" {
    #include "ilm_client.h"
    #include "ilm_control.h"
}

class TestBase
{
public:
    TestBase();
    virtual ~TestBase();
    t_ilm_surface getSurface();
    t_ilm_layer   getLayer();
protected:
    std::vector<wl_surface *> wlSurfaces;
    wl_display*    wlDisplay;

private:
    wl_registry*   wlRegistry;
    wl_compositor* wlCompositor;
    t_ilm_surface  nextSurfaceId;
    t_ilm_layer    nextLayerId;
};
