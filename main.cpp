#include <iostream>
#include <armadillo>
#include <gtk/gtk.h>

#include "bkmaxflow.hpp"

using namespace arma;
using namespace std;

namespace my_app
{
    GtkWidget *window, *draw_area_1, *draw_area_2;
    GdkPixbuf *image_pix, *result_pix;

    bk_maxflow::graph_cut * cutter;

    struct image_
    {
        int N,M;
        cube rgb;
    } image;

    cube seeds_mat;
    mat results_mat;

    int rad = 10;

    static gboolean click_load(gpointer ud)
    {
        GtkWidget * dialog = gtk_file_chooser_dialog_new("Open File",GTK_WINDOW(window),
        GTK_FILE_CHOOSER_ACTION_OPEN,
        "Cancel",GTK_RESPONSE_CANCEL,
        "Open" , GTK_RESPONSE_ACCEPT, NULL);
        if(gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
        {
            char * filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

            image_pix = gdk_pixbuf_new_from_file(filename,NULL);
            int W,H,row_stride,n_chan;

            guchar * pixels;

            n_chan = gdk_pixbuf_get_n_channels(image_pix);

            W = gdk_pixbuf_get_width(image_pix);
            H = gdk_pixbuf_get_height(image_pix);
            image.N = H;
            image.M = W;

            row_stride = gdk_pixbuf_get_rowstride(image_pix);
            pixels = gdk_pixbuf_get_pixels(image_pix);

            if(n_chan >= 3)
            {
                image.rgb.resize(H,W,3);
                seeds_mat.resize(H,W,2);
                seeds_mat.fill(0.0);
            }

            for(int y=0;y<H;y++)
            for(int x=0;x<W;x++)
            {
                image.rgb(y,x,0) = pixels[y*row_stride + x*n_chan];
                image.rgb(y,x,1) = pixels[y*row_stride + x*n_chan + 1];
                image.rgb(y,x,2) = pixels[y*row_stride + x*n_chan + 2];
            }

            gtk_widget_show_all(window);
            gtk_widget_set_size_request(draw_area_1,W,H);
            gtk_widget_set_size_request(draw_area_2,0,0);
            gtk_widget_set_size_request(window,MIN(W,H),H);

            gtk_widget_destroy(dialog);
            return TRUE;
        }
        return FALSE;
    }

    static gboolean click_graph(gpointer ud)
    {
        int wind = 2;
        int v_size = (2*wind+1)*(2*wind+1)-1;
        cutter = new bk_maxflow::graph_cut();

        g_print("Adding nodes to graph...");
        cutter->add_nodes(image.N*image.M);
        g_print("Done \n");

        g_print("Adding edges to graph...");
        cutter->add_edges(v_size * image.N * image.M);
        g_print("Done \n");
        vec hist_red,hist_blue;
        ivec scale; scale << 8 << 8 << 8;
        ivec size = 256 / scale;

        hist_red = zeros<vec>(size(0)*size(1)*size(2));
        hist_blue = zeros<vec>(size(0)*size(1)*size(2));

        g_print("Preparing Histograms...");
        for(int u=0;u<image.N;u++)
            for(int v=0;v<image.M;v++)
            {
                int idx =   size(1) * size(2) * (int) (image.rgb(u,v,0) / scale(0) )
                          + size(2) * (int) (image.rgb(u,v,1) / scale(1) )
                          + (int) (image.rgb(u,v,2) / scale(2));
                if(seeds_mat(u,v,0) > 0)
                    hist_red(idx) += 1.0;
                if(seeds_mat(u,v,1) > 0)
                    hist_blue(idx) += 1.0;
            }

        hist_red /= sum(hist_red);
        hist_red = -log(hist_red);
        hist_blue /= sum(hist_blue);
        hist_blue = -log(hist_blue);
        g_print("Done\n");

        g_print("Finding variance...");
        double beta = 2.0;
        mat sigma = zeros<mat>(3,3);
        vec mu = zeros<vec>(3);
        double inv_NM = 1.0 / (image.N * image.M);

        for(int u=0;u<image.N;u++)
            for(int v=0;v<image.M;v++)
                mu += inv_NM*image.rgb.tube(u,v);
        cout<<mu<<endl;

        for(int u=0;u<image.N;u++)
            for(int v=0;v<image.M;v++)
            {
                vec aux;
                aux = vec(image.rgb.tube(u,v)) - mu;
                sigma += inv_NM * aux * aux.t();
            }
        cout<<sigma<<endl;
        sigma = inv(sigma);

        g_print("Done\n");

        g_print("Preparing weights...");
        double K = 1.0;
        for(int u=0;u<image.N;u++)
            for(int v=0;v<image.M;v++)
            {
                double sum = 1.0;
                for(int i=-wind;i<=wind;i++)
                    for(int j=-wind;j<=wind;j++)
                        if(i != 0 || j != 0)
                        {
                            if(u + i < 0 || v + j < 0) continue;
                            if(u + i >= image.N || v + j >= image.M) continue;

                            int p_idx = image.M*u + v;
                            int q_idx = image.M*(u+i)+ v+j;

                            vec aux = (image.rgb.tube(u,v) - image.rgb.tube(u+i,v+j));

                            double w_uv = exp(-beta*as_scalar(aux.t()*sigma*aux)) / sqrt(i*i + j*j);
                            sum += w_uv;

                            if(!cutter->connected(p_idx,q_idx))
                                cutter->set_edge_weight(p_idx,q_idx,w_uv,w_uv);
                        }
                K = K < sum ? sum : K;
            }

        double lambda = 1e5;
        for(int u=0;u<image.N;u++)
            for(int v=0;v<image.M;v++)
            {
                int p_idx = image.M * u + v;
                if(seeds_mat(u,v,0) > 0)
                    cutter->set_terminal_weights(p_idx,K,0);
                else if(seeds_mat(u,v,1) > 0)
                    cutter->set_terminal_weights(p_idx,0,K);
                else
                {
                    int idx =   size(1) * size(2) * (int) (image.rgb(u,v,0) / scale(0) )
                              + size(2) * (int) (image.rgb(u,v,1) / scale(1) )
                              + (int) (image.rgb(u,v,2) / scale(2));
                    cutter->set_terminal_weights(p_idx,lambda*hist_blue(idx),lambda*hist_red(idx));
                }
            }
        g_print("Done\n");

        g_print("Doing BK Maxflow...");
        results_mat.resize(image.N,image.M);
        std::vector<int> new_nodes;
        double flow = cutter->get_maxflow(false,new_nodes);
        g_print("Done\n");

        for(int u=0;u<image.N;u++)
            for(int v=0;v<image.M;v++)
            {
                int p_idx = image.M*u + v;
                if(cutter->get_terminal(p_idx) == bk_maxflow::FOREGROUND)
                {
                    results_mat(u,v) = 1;
                }
                else results_mat(u,v) = 0;
            }

        result_pix = gdk_pixbuf_copy(image_pix);
        guchar * pix = gdk_pixbuf_get_pixels(result_pix);

        int row_stride = gdk_pixbuf_get_rowstride(result_pix);
        int n_chan = gdk_pixbuf_get_n_channels(result_pix);

        for(int u=0;u<image.N;u++)
            for(int v=0;v<image.M;v++)
            {
            pix[u*row_stride + v*n_chan] = 0;
            pix[u*row_stride + v*n_chan + 1] = 0;
            pix[u*row_stride + v*n_chan + 2] = 0;
            if(results_mat(u,v)>0)
                {
                pix[u*row_stride + v*n_chan] = 255;
                pix[u*row_stride + v*n_chan + 1] = 255;
                pix[u*row_stride + v*n_chan + 2] = 255;
                }
            }

        gtk_widget_set_size_request(draw_area_2,image.M,image.N);
        delete cutter;
        return TRUE;
    }

    static gboolean click_over_image(GtkWidget * widget,GdkEventButton * event,gpointer data)
    {
        g_print("click: %f %f ",event->x,event->y);
        if(event->button == 1)
            {
            g_print("left\n");
            }
        if(event->button == 3)
            {
            g_print("right\n");
            }
        else
            g_print("\n");
        return TRUE;
    }

    gboolean mouse_over_image(GtkWidget *widget,GdkEventMotion *event,gpointer data)
    {
        GdkModifierType state;
        GdkDeviceManager * device_manager;
        GdkDevice * pointer;
        device_manager = gdk_display_get_device_manager(gtk_widget_get_display(widget));
        pointer = gdk_device_manager_get_client_pointer(device_manager);

        int x, y;

        gdk_window_get_device_position (event->window, pointer,&x, &y, &state);

        if(draw_area_1 == NULL) return FALSE;
        if(image_pix == NULL) return FALSE;
        guchar * pix = gdk_pixbuf_get_pixels(image_pix);

        int row_stride = gdk_pixbuf_get_rowstride(image_pix);
        int n_chan = gdk_pixbuf_get_n_channels(image_pix);

        if(state & GDK_BUTTON1_MASK)
        {
            int inf_x = MAX(0,x-rad);
            int inf_y = MAX(0,y-rad);
            int sup_x = MIN(image.M-1,x+rad);
            int sup_y = MIN(image.N-1,y+rad);
            for(int u = inf_x;u <= sup_x;u++)
                for(int v = inf_y;v <= sup_y;v++)
                    if((u-x)*(u-x) + (v-y)*(v-y) <= rad*rad)
                    {
                        seeds_mat(v,u,0) = 255;
                        seeds_mat(v,u,1) = 0;
                        pix[v*row_stride + u*n_chan] =  255;
                    }
        }
        else if(state & GDK_BUTTON3_MASK)
        {
            int inf_x = MAX(0,x-rad);
            int inf_y = MAX(0,y-rad);
            int sup_x = MIN(image.M-1,x+rad);
            int sup_y = MIN(image.N-1,y+rad);
            for(int u = inf_x;u <= sup_x;u++)
                for(int v = inf_y;v <= sup_y;v++)
                    if((u-x)*(u-x) + (v-y)*(v-y) <= rad*rad)
                    {
                        seeds_mat(v,u,0) = 0;
                        seeds_mat(v,u,1) = 255;
                        pix[v*row_stride + u*n_chan+2] =  255;
                    }
        }
        gtk_widget_queue_draw(draw_area_1);
        return TRUE;
    }

    static gboolean da_expose1(GtkWidget *widget, GdkEventExpose *event)
    {
        if(image_pix)
        {
            cairo_t * cr = gdk_cairo_create(gtk_widget_get_window(widget));
            gdk_cairo_set_source_pixbuf(cr,image_pix,0,0);
            cairo_paint(cr);
            cairo_fill(cr);
            cairo_destroy(cr);
            gtk_widget_queue_draw(widget);
            return TRUE;
        }
        return FALSE;
    }

    static gboolean da_expose2(GtkWidget *widget, GdkEventExpose *event)
    {
        if(result_pix)
        {
            cairo_t * cr = gdk_cairo_create(gtk_widget_get_window(widget));
            gdk_cairo_set_source_pixbuf(cr,result_pix,0,0);
            cairo_paint(cr);
            cairo_fill(cr);
            cairo_destroy(cr);
            gtk_widget_queue_draw(widget);
            return TRUE;
        }
        return FALSE;
    }

    void init()
    {
        window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

        gtk_window_set_title(GTK_WINDOW(window),"Graphcut");
        g_signal_connect(G_OBJECT(window),"destroy",G_CALLBACK(gtk_main_quit),0);
        gtk_window_set_resizable(GTK_WINDOW(window),false);

        GtkWidget * grid = gtk_grid_new();

            GtkWidget * button = gtk_button_new_with_label("Load");
            g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(click_load),0);
            gtk_grid_attach(GTK_GRID(grid),button,0,0,1,1);

            button = gtk_button_new_with_label("Graphcut");
            g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(click_graph),0);
            gtk_grid_attach(GTK_GRID(grid),button,0,1,1,1);

            button = gtk_event_box_new();
            gtk_widget_set_size_request(button,1,1);
            g_signal_connect(button,"button-press-event",G_CALLBACK(click_over_image),NULL);
            g_signal_connect(button,"motion-notify-event",G_CALLBACK(mouse_over_image),NULL);

            draw_area_1 = gtk_drawing_area_new();
            gtk_widget_set_size_request(draw_area_1,500,500);
            g_signal_connect(draw_area_1,"draw",G_CALLBACK(da_expose1),NULL);
            gtk_container_add(GTK_CONTAINER(button),draw_area_1);

            draw_area_2 = gtk_drawing_area_new();
            gtk_widget_set_size_request(draw_area_2,0,0);
            g_signal_connect(draw_area_2,"draw",G_CALLBACK(da_expose2),NULL);

        gtk_grid_attach(GTK_GRID(grid),button ,3,0,2,2);
        gtk_grid_attach(GTK_GRID(grid),draw_area_2,5,0,2,2);
        gtk_container_add(GTK_CONTAINER(window),grid);

        image.rgb.resize(0,0,0);
        image.M = image.N = 0;
        gtk_widget_show_all(window);
    }

};

int main(int n_params,char ** params)
{
    gtk_init(&n_params,&params);

    my_app::init();

    gtk_main();
    return 0;
}
