
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
protected:
    std::vector<wl_surface *> wlSurfaces;
    wl_display*    wlDisplay;

private:
    wl_registry*   wlRegistry;
    wl_compositor* wlCompositor;
    static t_ilm_surface  nextSurfaceId;
    static t_ilm_surface  maxSurfaceId;
    static t_ilm_surface  startSurfaceId;
    static t_ilm_layer    nextLayerId;
    static t_ilm_layer    maxLayerId;
    static t_ilm_layer    startLayerId;

public:
    static t_ilm_surface getSurface()
    {
        if (nextSurfaceId >= (startSurfaceId + maxSurfaceId))
        {
            nextSurfaceId = startSurfaceId;
        }
        else
        {
            nextSurfaceId++;
        }

        return nextSurfaceId;
    }
    static t_ilm_layer getLayer()
    {
        if (nextLayerId >= (startLayerId + maxLayerId))
        {
            nextLayerId = startLayerId;
        }
        else
        {
            nextLayerId++;
        }

        return nextLayerId;
    }
    static t_ilm_layer getStartLayerId() { return startLayerId; }
    static t_ilm_layer getEndLayerId() { return (startLayerId + maxLayerId); }
    static t_ilm_layer getStartSurfaceId() { return startSurfaceId; }
    static t_ilm_layer getEndSurfaceId() { return (startSurfaceId + maxSurfaceId); }
    static void setStartSurfaceId(t_ilm_surface surfaceID) { startSurfaceId = surfaceID; nextSurfaceId = surfaceID;}
    static void setStartLayerId(t_ilm_layer layerID) { startLayerId = layerID; nextLayerId = layerID;}
    static void setMaxSurfaceIds(t_ilm_surface maxSurfaceID) { maxSurfaceId = maxSurfaceID; }
    static void setMaxLayerIds(t_ilm_layer maxLayerID) { maxLayerId = maxLayerID; }
};
