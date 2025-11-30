use cgmath::prelude::*;
use three_d::{degrees, Camera, Event, MouseButton};


pub fn orbit_control_update_camera(
    camera: &mut Camera, events: &mut [Event], min_dist: f32, max_dist: f32
)
    -> bool
{
    let target = camera.target();
    let mut change = false;

    for event in events.iter_mut() {
        match event {
            Event::MouseMotion { delta, button, handled, .. } => {
                if let Some(button) = button {
                    if !*handled {
                        let (x, y) = *delta;

                        match button {
                            MouseButton::Left => {
                                camera.rotate_around(target, 0.1 * x, 0.1 * y);
                            }
                            MouseButton::Middle => {
                                let cam_up = camera.up_orthogonal();
                                let cam_right = camera.right_direction();
                                camera.translate(0.01 * (y * cam_up - x * cam_right));
                            }
                            MouseButton::Right => {
                                camera.roll(degrees(0.1 * x));
                            }
                        }
                        *handled = true;
                        change = true;
                    }
                }
            }
            Event::MouseWheel { delta, handled, .. } => {
                if !*handled {
                    let dist = target.distance(camera.position());
                    camera.zoom_towards(target, 0.001 * delta.1 * dist, min_dist, max_dist);
                    *handled = true;
                    change = true;
                }
            }
            _ => {}
        }
    }
    change
}
