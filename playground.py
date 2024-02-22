import subprocess
import tensorflow as tf
options = tf.profiler.experimental.ProfilerOptions(host_tracer_level = 3,
                                                   python_tracer_level = 1,
                                                   device_tracer_level = 1)

logdir = "./examples/models/logs/"

tf.profiler.experimental.start(logdir, options = options)
print('fick dich')
tf.profiler.experimental.stop()

subprocess.run(['tensorboard', '--logdir', logdir])
