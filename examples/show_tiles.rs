use cgmath::prelude::*;
use cgmath::{point3, vec3, vec4, Point3};
use rust_dsymbols::dsets::DSet;
use rust_dsymbols::geometry::vec_matrix::VecMatrix;
use three_d::Mat4;

use rust_dsymbols::delaney3d::pseudo_toroidal_cover;
use rust_dsymbols::display::mesh::{ItemType, Mesh, decompose_mesh, scaled_mesh};
use rust_dsymbols::dsyms::PartialDSym;
use rust_dsymbols::tilings::{Skeleton, chamber_positions, gram_matrix, invariant_basis, tile_surfaces};


struct Options {
    tile_scale: f64,
    edge_radius: f64,
    vertex_color: three_d::Srgba,
    edge_color: three_d::Srgba,
    face_color: three_d::Srgba,
}


impl Default for Options {
    fn default() -> Self {
        Self {
            tile_scale: 0.75,
            edge_radius: 0.05,
            vertex_color: three_d::Srgba::BLACK,
            edge_color: three_d::Srgba::BLUE,
            face_color: three_d::Srgba::RED,
        }
    }
}


fn main() {
    #[cfg(feature = "pprof")]
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build().unwrap();

    // On the web, this creates a canvas instead.
    let window = three_d::Window::new(three_d::WindowSettings {
        title: "Rust 3d Test".to_string(),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();

    let mut gui = three_d::GUI::new(&context);

    let mut camera = three_d::Camera::new_perspective(
        window.viewport(),
        vec3(0.0, 2.0, 8.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        three_d::degrees(25.0),
        0.1,
        1000.0,
    );

    let control = rust_dsymbols::display::controls::OrbitControl::new(1.0, 1000.0);

    //let ds_spec = "<1.1:2 3:1 2,1 2,1 2,2:3 3,4 3,4>";
    let ds_spec = "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>";

    let options = Default::default();

    let models = build_models(&context, ds_spec, &options);

    let sun_dir = vec4(1.0, -1.0, -1.0, 0.0);

    let mut sun = three_d::DirectionalLight::new(
        &context,
        2.0,
        three_d::Srgba::WHITE,
        sun_dir.truncate()
    );

    let ambient = three_d::AmbientLight::new(
        &context,
        0.1,
        three_d::Srgba::WHITE,
    );

    let mut shadows_enabled = false;

    #[cfg(feature = "pprof")]
    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };

    window.render_loop(move |mut frame_input| {
        let mut panel_width = 0.0;
        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |gui_context| {
                use three_d::egui::*;
                SidePanel::left("side_panel").show(gui_context, |ui| {
                    ui.heading("Settings");

                    ui.label("Appearance");

                    if ui.checkbox(&mut shadows_enabled, "Shadows On").clicked() &&
                        !shadows_enabled
                    {
                        sun.clear_shadow_map();
                    }
                });
                panel_width = gui_context.used_rect().width();
            },
        );

        let viewport = three_d::Viewport {
            x: (panel_width * frame_input.device_pixel_ratio) as i32,
            y: 0,
            width: frame_input.viewport.width
                - (panel_width * frame_input.device_pixel_ratio) as u32,
            height: frame_input.viewport.height,
        };
        camera.set_viewport(viewport);
        control.handle_events(&mut camera, &mut frame_input.events);

        // Moves the sun around the objects in sync with the camera.
        sun.direction = (camera.view().invert().unwrap() * sun_dir).truncate();

        if shadows_enabled {
            sun.generate_shadow_map(2048, &models);
        }

        let screen = frame_input.screen();

        screen.clear(three_d::ClearState::color_and_depth(0.8, 0.8, 0.8, 1.0, 1.0));
        screen.render(&camera, &models, &[&sun, &ambient]);
        screen.write(|| gui.render()).unwrap();

        // Ensures a valid return value.
        three_d::FrameOutput::default()
    });
}


fn build_models(context: &three_d::Context, ds_spec: &str, options: &Options)
    -> Vec<three_d::Gm<three_d::InstancedMesh, three_d::PhysicalMaterial>>
{
    let instances = three_d::Instances {
        transformations: vec![
            Mat4::from_scale(1.0),
        ],
        ..Default::default()
    };

    tiles(ds_spec).iter().flat_map(|tile_mesh|
        decompose_mesh(
                &scaled_mesh(tile_mesh, options.tile_scale),
                options.edge_radius
            )
            .iter()
            .map(|(part_mesh, item_type)| {
                let color = match item_type {
                    ItemType::Vertex => options.vertex_color,
                    ItemType::Edge => options.edge_color,
                    ItemType::Face => options.face_color,
                };

                three_d::Gm::new(
                    three_d::InstancedMesh::new(
                        context,
                        &instances,
                        &part_mesh.to_cpu_mesh()
                    ),
                    three_d::PhysicalMaterial {
                        albedo: color,
                        metallic: 0.0,
                        roughness: 0.5,
                        ..Default::default()
                    }
                )
            })
            .collect::<Vec<_>>()
        )
        .collect()
}


fn tiles(spec: &str) -> Vec<Mesh<Point3<f64>>> {
    let ds = spec.parse::<PartialDSym>().unwrap();
    let cov = pseudo_toroidal_cover(&ds).unwrap();
    let skel = Skeleton::of(&cov);
    let basis = invariant_basis(&gram_matrix(&ds, &cov, &skel).unwrap()).transpose();
    let pos = skel.graph.vertices().iter()
        .map(|&v| (v, skel.graph.position(v)))
        .collect();

    let (scale, shift) = scale_and_shift(&cov, &skel, &basis);
    let reps = cov.orbit_reps([0, 1, 2], cov.elements());

    tile_surfaces(&cov, &skel, &pos, reps).iter().map(|(vertices, faces)| {
        let vs = vertices.iter().map(|v| {
            let v = scale * &basis * v.to_f64().unwrap() + &shift;
            point3(v[(0, 0)], v[(1, 0)], v[(2, 0)])
        });

        Mesh::from_oriented_faces(vs, faces.clone()).unwrap()
    }).collect()
}


fn scale_and_shift(cov: &PartialDSym, skel: &Skeleton, basis: &VecMatrix<f64>)
    -> (f64, VecMatrix<f64>)
{
    let ch_pos = chamber_positions(cov, skel);
    let mut s = 0.0;
    let mut n = 0;

    for d in  cov.orbit_reps([0, 2, 3], cov.elements()) {
        let v = 2.0 * basis * (&ch_pos[&(d, 1)] - &ch_pos[&(d, 0)]).to_f64().unwrap();
        s += v.norm();
        n += 1;
    };

    let scale = n as f64 / s;
    let shift = -scale * basis * ch_pos[&(1, 3)].to_f64().unwrap();

    (scale, shift)
}


trait ToCpuMesh {
    fn to_cpu_mesh(&self) -> three_d::CpuMesh;
}


impl ToCpuMesh for Mesh<Point3<f64>> {
    fn to_cpu_mesh(&self) -> three_d::CpuMesh {
        let trimesh = self.triangulate().unwrap();

        let positions: Vec<_> = trimesh.vertices().iter()
            .map(|&p| p.to_vec())
            .collect();

        let indices: Vec<_> = trimesh.face_indices().iter()
            .flatten()
            .map(|&x| x as u32)
            .collect();

        let mut mesh = three_d::CpuMesh {
            positions: three_d::Positions::F64(positions),
            indices: three_d::Indices::U32(indices),
            ..Default::default()
        };
        mesh.compute_normals();
        mesh
    }
}
