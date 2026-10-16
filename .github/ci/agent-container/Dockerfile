# The default image of the agent's container backend
# (agent.bash_isolation = "container", delfin/agent/container_exec.py
# DEFAULT_OCI_IMAGE). Nothing builds this file. It exists so Dependabot's
# docker ecosystem, which reads Dockerfiles and not Python, proposes the new
# digest when python:3.11-slim moves; such a pull request fails
# tests/test_what_the_agent_runs_on_is_named_exactly.py until
# DEFAULT_OCI_IMAGE carries the same digest, so the image the agent runs is
# changed by a person who read the update.
#
# Digest: the OCI image index of docker.io/library/python:3.11-slim, read
# from the registry on 2026-10-16 (Docker-Content-Digest of the manifest,
# equal to the sha256 of its body).
FROM python:3.11-slim@sha256:e88e9763f943ec1834f992a4b51e0f24500486803e8bc534e5767af9ea65f6ce
