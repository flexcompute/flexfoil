# Running the FlexFoil web UI (remote)

The dev server is Vite (port **5173**, binds `localhost`). The multi-element
high-lift design tool is at `/highlift-preview.html`; the main app is at `/`.

## 1. Launch on the remote machine (`014-v100-dev`)

```bash
export PATH="$HOME/.nvm/versions/node/v22.22.3/bin:$PATH"   # Node 22 (Vite needs ≥18)
cd /home/qiqi/flexcompute/flexfoil/flexfoil-ui
npm install        # first time only — node_modules is already present
npm run dev        # serves http://localhost:5173  (leave running)
```

## 2. Open an SSH tunnel from your **local** machine

```bash
ssh -N -L 5173:localhost:5173 qiqi@014-v100-dev
```

- `-N` = just forward, run no remote command. Keep this open while viewing.
- If you reach the box through a jump host or a different name/IP, use that ssh
  target — only the `-L 5173:localhost:5173` forward matters.
- Local port 5173 busy? Use e.g. `-L 8080:localhost:5173` and browse `:8080`.

## 3. Open in your local browser

- High-lift multi-element tool → http://localhost:5173/highlift-preview.html
- Main FlexFoil app → http://localhost:5173/

Hot reload (HMR) works through the tunnel — edits under `flexfoil-ui/src/` refresh
live.
